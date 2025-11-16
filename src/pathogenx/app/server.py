import pandas as pd
from scipy.sparse import coo_matrix
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
from pathogenx.io import GenotypeFile, MetaFile, DistFile, _GENOTYPE_FLAVOURS, _META_FLAVOURS, _DIST_FLAVOURS
from pathogenx.dataset import Dataset
from shinywidgets import render_widget
from pathogenx.calculators import PrevalenceCalculator, PrevalenceResult

from pathogenx.app import DATA, RESMAP, X_AXIS_VARIABLES
from pathogenx.app.plotters import (PrevalencePlotter, StrataPlotter, SummaryBarPlotter, CoveragePlotter, MapPlotter,
                                 merge_prevalence_figs)


def upload_server(input: Inputs, output: Outputs, session: Session):
    # Reactive container for user uploaded files to be loaded ----------------------
    @reactive.Calc
    def reactive_genotypes() -> pd.DataFrame | None:
        files = input.genotype_()
        if not files:
            return None
        file_info = [{"name": f["name"], "size": f["size"], "type": f["type"]} for f in files]
        return pd.DataFrame(file_info)

    @reactive.Calc
    def reactive_metadata() -> pd.DataFrame | None:
        files = input.upload()
        if not files:
            return None
        file_info = [{"name": f["name"], "size": f["size"], "type": f["type"]} for f in files]
        return pd.DataFrame(file_info)

    @reactive.Calc
    def reactive_distances() -> DistFile | None:
        files = input.upload()
        if not files:
            return None

        file_info = [{"name": f["name"], "size": f["size"], "type": f["type"]} for f in files]
        return pd.DataFrame(file_info)


    # Loads all data uploaded by user and calculates clusters ----------------------
    @reactive.Calc
    def reactive_dataset() -> Dataset | None:
        files_df = reactive_files()
        if files_df is None or files_df.empty:
            return None

        n_collections = 1 # Placeholder for len(files_df["collection"].unique())
        ui.notification_show(f"Loading {n_collections} collection(s)...", type="message")



    # User data to be filtered and used for prevalences ----------------------------
    @reactive.Calc
    def reactive_data() -> pd.DataFrame | None:
        d = reactive_loaded_data()
        if d is None or d.empty:
            return None

        # Wait for inputs to be available
        req_inputs = [
            input.study_selector(),
            input.year_selector(),
            input.region_selector(),
        ]
        if not all(i is not None for i in req_inputs):
            return None

        # Add spatio-temporal data if selected
        # In this Python version, we assume columns are already named correctly
        # or would be renamed in the read_pw function.

        # Apply filters
        year_min, year_max = input.year_selector()

        amr_map = {"ESBL+": 1, "Carbapenemase+": 2}
        min_resistance = amr_map.get(input.amr_selector(), 0)

        filtered_d = d[
            d["Year"].between(year_min, year_max) &
            d["Study"].isin(input.study_selector()) &
            d["Region"].isin(input.region_selector()) &
            (d["resistance_score"] >= min_resistance)
        ]

        if filtered_d.empty:
            ui.notification_show("No strains matching filter", type="warning")

        return filtered_d

    # Reactive variables -----------------------------------------------------------
    @reactive.Calc
    def reactive_metadata_variables() -> list[str]:
        d = reactive_loaded_data()
        if d is None or d.empty:
            return []

        known_cols = ["id", "K_locus", "O_type", "Study", "Cluster", "Cluster_size", "resistance_score", "ST"] # Example
        return [col for col in d.columns if col not in known_cols]

    @reactive.Effect
    def _update_geo_pickers():
        meta_vars = reactive_metadata_variables()
        choices = meta_vars if meta_vars else []
        ui.update_select(session, "year_variable", choices=choices)
        ui.update_select(session, "region_variable", choices=choices)
        ui.update_select(session, "country_variable", choices=choices)

    # Update study selector
    @reactive.Effect
    @reactive.event(reactive_loaded_data)
    def _update_study_selector():
        d = reactive_loaded_data()
        if d is None or d.empty:
            studies = []
        else:
            studies = sorted(d["Study"].unique().tolist())

        ui.update_select(session, "study_selector", choices=studies, selected=studies)

    # Update region selector
    @reactive.Effect
    @reactive.event(reactive_loaded_data)
    def _update_region_selector():
        d = reactive_loaded_data()
        if d is None or d.empty or "Region" not in d.columns:
            regions = []
        else:
            regions = sorted(d["Region"].dropna().unique().tolist())

        ui.update_select(session, "region_selector", choices=regions, selected=regions)

    # Update year slider
    @reactive.Effect
    @reactive.event(reactive_loaded_data)
    def _update_year_slider():
        d = reactive_loaded_data()
        if d is None or d.empty or "Year" not in d.columns:
            min_y, max_y = 1900, 2100
        else:
            years = d["Year"].dropna()
            min_y, max_y = int(years.min()), int(years.max())

        ui.update_slider("year_selector", min=min_y, max=max_y, value=(min_y, max_y))

    # Data selector reset ----------------------------------------------------------
    @reactive.Effect
    @reactive.event(input.data_reset)
    def _reset_data_filters():
        ui.notification_show("Resetting data filters", type="message")
        ui.update_select(session, "amr_selector", selected="All samples")
        ui.update_select(session, "year_variable", selected=[])
        ui.update_select(session, "region_variable", selected=[])
        ui.update_select(session, "country_variable", selected=[])

        # Re-trigger updates for study, region, year
        _update_study_selector()
        _update_region_selector()
        _update_year_slider()


    # Control panels ---------------------------------------------------------------
    @reactive.Effect
    @reactive.event(reactive_loaded_data)
    def _toggle_panels_on_load():
        d = reactive_loaded_data()

        # Always remove panels first to handle re-uploads
        ui.remove_ui(selector="#prevalence_panel", immediate=True)
        ui.remove_ui(selector="#coverage_panel", immediate=True)
        ui.remove_ui(selector="#dataframe_panel", immediate=True)

        if d is None or d.empty:
            ui.toggle_sidebar("sidebar", open=False)
            ui.accordion_panel_open("accordion", "upload_panel")
        else:
            ui.accordion_panel_close("accordion", "upload_panel")

            # Insert panels
            ui.insert_ui(
                selector="div#accordion",
                where="beforeEnd",
                ui=prevalence_panel,
            )
            ui.insert_ui(
                selector="div#accordion",
                where="beforeEnd",
                ui=coverage_panel,
            )
            ui.insert_ui(
                selector="div#accordion",
                where="beforeEnd",
                ui=dataframe_panel,
            )
            ui.toggle_sidebar("sidebar", open=True)

    # Reactive antigen lists -------------------------------------------------------
    @reactive.Effect
    @reactive.event(reactive_loaded_data, input.antigen_selector, input.vaccine_valency)
    def _update_antigen_order_inputs():
        d = reactive_loaded_data()
        if d is None or d.empty:
            return

        antigen_col = input.antigen_selector()
        if antigen_col not in d.columns:
            return

        # This should be based on prevalence, but for now, just use value_counts
        all_antigens = d[antigen_col].value_counts().index.tolist()

        valency = min(input.vaccine_valency(), len(all_antigens))
        if valency == 0:
            return

        default_antigens = all_antigens[:valency]
        available_antigens = all_antigens[valency:]


    # Output dataframe -------------------------------------------------------------
    @render.data_frame
    def dataframe():
        df = reactive_data()
        if df is None or df.empty:
            return
        return render.DataGrid(df.reset_index(drop=True))

    # Output summary ---------------------------------------------------------------
    @render.text
    def summary():
        loaded = reactive_loaded_data()
        if loaded is None or loaded.empty:
            return ""

        filtered = reactive_data()
        if filtered is None or filtered.empty:
            return f"Samples: 0/{len(loaded)}"

        # This is a simplified version of your R summary string
        return (
            f"Samples: {len(filtered)}/{len(loaded)}; "
            f"Studies: {len(input.study_selector() or [])}/{len(loaded['Study'].unique())}; "
            f"Regions: {len(input.region_selector() or [])}/{len(loaded['Region'].unique())}; "
            f"Years: {input.year_selector()[0]}-{input.year_selector()[1]}; "
            f"Resistance: {input.amr_selector()}"
        )

def server(input: Inputs, output:Outputs, session: Session):
    @reactive.Calc
    def reactive_data():
        """Filters the main dataframe based on UI inputs."""
        filtered_df = DATA[
            (DATA['Year'].between(input.year_selector()[0], input.year_selector()[1])) &
            (DATA['Study'].isin(input.study_selector())) &
            (DATA['Region'].isin(input.region_selector())) &
            (DATA['resistance_score'] >= RESMAP.get(input.amr_selector(), 0)) &
            (DATA['Mortality'] == 1 if input.outcome_selector() != 'All samples' else True)
            ]
        if filtered_df.empty:
            ui.notification_show("No strains matching filter", type="warning")
        return filtered_df

    @reactive.Calc
    def reactive_prevalence() -> PrevalenceResult:
        """Calculates overall prevalence for the selected genotype."""
        return PrevalenceCalculator(
            stratify_by=[input.genotype_selector()],
            adjust_for=['Cluster'],
            n_distinct=[v for v in X_AXIS_VARIABLES if v not in {
                input.genotype_selector(), input.heatmap_x()}]
        ).calculate(reactive_data())

    @reactive.Calc
    def reactive_prevalence_stratified() -> PrevalenceResult:
        """Calculates prevalence stratified by a second variable for the heatmap."""
        return PrevalenceCalculator(
            stratify_by=[input.genotype_selector(), input.heatmap_x()],
            adjust_for=['Cluster'],
            denominator=(
                input.genotype_selector()
                if input.heatmap_swap_denominator()
                else input.heatmap_x()
            )
        ).calculate(reactive_data())

    @reactive.Calc
    def reactive_prevalence_coverage() -> PrevalenceResult:
        """Calculates prevalence stratified by a second variable for the heatmap."""
        return PrevalenceCalculator(
            stratify_by=[input.coverage_group(), input.genotype_selector()],
            adjust_for=['Cluster']
        ).calculate(reactive_data())

    @reactive.Calc
    def genotype_list() -> list[str]:
        """Determines the list of genotypes to display based on user selections."""
        if len(result := reactive_prevalence()) == 0:
            return []
        valency = min(input.vaccine_valency(), len(result))
        return result.data[input.genotype_selector()].iloc[:valency].tolist()

    @output
    @render_widget
    def global():
        """Renders the main combined plot (pyramid, heatmap, bars)."""
        # Get reactive variables
        r1 = reactive_prevalence()
        r2 = reactive_prevalence_stratified()
        y_order = input.genotype_selected()
        if len(r1) == 0 or len(r2) == 0 or not y_order:
            return None
        # Init plotters
        p1 = PrevalencePlotter(y_order=y_order)
        p2 = StrataPlotter(y_order=y_order, max_x=input.heatmap_num_x())
        p3 = SummaryBarPlotter(y_order=y_order, fill_by=input.bars_x())
        # Render plots
        return merge_prevalence_figs(p1.plot(r1), p2.plot(r2), p3.plot(r1))

    @output
    @render_widget
    def coverage_plot():
        """Renders the main combined plot (pyramid, heatmap, bars)."""
        if len(r := reactive_prevalence_coverage()) == 0 or not (x_order := genotype_list()):
            return None
        return CoveragePlotter(x_order=x_order).plot(r)

    @output
    @render_widget
    def map():
        """Renders the main combined plot (pyramid, heatmap, bars)."""
        if len(r := reactive_prevalence_coverage()) == 0:
            return None
        return  MapPlotter().plot(r)

    @output
    @render.text
    def summary():
        """
        Returns a summary of the filtered data
        """
        return (f"Samples: {len(reactive_data())}/{len(DATA)}; "
            f"Studies: {len(input.study_selector())}/{len(DATA['Study'].unique())}; "
            f"Regions: {len(input.region_selector())}/{len(DATA['Region'].unique())}; "
            f"Years: {input.year_selector()[1] - input.year_selector()[0] + 1}/"
            f"{DATA['Year'].max() - DATA['Year'].min() + 1}; "
            f"Resistance: {input.amr_selector()}; "
            f"Outcome: {input.outcome_selector()}")
