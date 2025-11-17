import pandas as pd
from scipy.sparse import coo_matrix
from pathlib import Path
from shiny import Inputs, Outputs, Session, reactive, render, ui
from pathogenx.io import GenotypeFile, MetaFile, DistFile, _GENOTYPE_FLAVOURS, _META_FLAVOURS, _DIST_FLAVOURS
from pathogenx.dataset import Dataset
from pathogenx.calculators import PrevalenceCalculator, PrevalenceResult
from pathogenx.app.plotters import (PrevalencePlotter, StrataPlotter, SummaryBarPlotter, CoveragePlotter, MapPlotter,
                                 merge_prevalence_figs)


def main_server(input: Inputs, output: Outputs, session: Session):
    # Reactive container for user uploaded files to be loaded ----------------------
    @reactive.calc
    def reactive_genotypes() -> pd.DataFrame | None:
        file_infos = input.genotype_upload()
        if not file_infos:
            return None
        f = file_infos[0]  # Get the first file info dictionary
        genotype_file = GenotypeFile.from_flavour(Path(f["datapath"]), input.genotype_flavour())
        try:
            data = genotype_file.load()
            ui.notification_show('Successfully loaded genotype file', type='message')
            return data
        except Exception as e:
            ui.notification_show(f"Error loading genotype file: {e}", type="error")
            return None

    @reactive.calc
    def reactive_metadata() -> pd.DataFrame | None:
        file_infos = input.metadata_upload()
        if not file_infos:
            return None
        f = file_infos[0]
        metadata_file = MetaFile.from_flavour(Path(f["datapath"]), input.metadata_flavour())
        try:
            data = metadata_file.load()
            ui.notification_show('Successfully loaded metadata file', type='message')
            return data
        except Exception as e:
            ui.notification_show(f"Error loading metadata file: {e}", type="error")
            return None

    @reactive.calc
    def reactive_distances() -> tuple[coo_matrix, list[str]] | None:
        file_infos = input.distance_upload()
        if not file_infos:
            return None
        f = file_infos[0]
        distance_file = DistFile.from_flavour(Path(f["datapath"]), input.distance_flavour())
        try:
            data = distance_file.load()
            ui.notification_show('Successfully loaded distance file', type='message')
            return data
        except Exception as e:
            ui.notification_show(f"Error loading distance file: {e}", type="error")
            return None

    # Create a reactive value to hold the dataset state
    reactive_dataset = reactive.Value[Dataset | None](None)

    @reactive.effect
    @reactive.event(reactive_genotypes, reactive_metadata, reactive_distances)
    def _create_dataset():
        """
        This effect runs whenever the underlying file data changes.
        It creates a new Dataset object and sets it as the reactive value.
        """
        if (genotypes := reactive_genotypes()) is None or genotypes.empty:
            reactive_dataset.set(None)
            return
        dataset = Dataset(genotypes, reactive_metadata(), reactive_distances())
        reactive_dataset.set(dataset)

    @reactive.event(input.calc_clusters)
    def _calculate_clusters():
        """When the button is clicked, this modifies the dataset in-place."""
        d = reactive_dataset.get()
        if d is not None and d.distances is not None:
            ui.notification_show('Calculating clusters...')
            d.calculate_clusters(input.cluster_method(), input.snp_distance())
            # Since we modified the object, we need to trigger its dependents to re-run.
            # We do this by "setting" the value to itself.
            reactive_dataset.set(d)
        else:
            ui.notification_show('Cannot calculate clusters without a distance file.', type='error')


    @reactive.event(input.reset_uploads)
    def _reset_uploads():
        ui.notification_show('Resetting uploads', type='message')


    # User data to be filtered and used for prevalences ----------------------------
    @reactive.calc
    def reactive_data() -> pd.DataFrame | None:
        if (d := reactive_dataset.get()) is None or len(d) == 0:
            return None

        # # Wait for inputs to be available
        # req_inputs = [
        #     input.study_selector(),
        #     input.year_selector(),
        #     input.region_selector(),
        # ]
        # if not all(i is not None for i in req_inputs):
        #     return None
        # 
        # # Add spatio-temporal data if selected
        # # In this Python version, we assume columns are already named correctly
        # # or would be renamed in the read_pw function.
        # 
        # # Apply filters
        # year_min, year_max = input.year_selector()
        # 
        # amr_map = {"ESBL+": 1, "Carbapenemase+": 2}
        # min_resistance = amr_map.get(input.amr_selector(), 0)
        # 
        # filtered_d = d[
        #     d["Year"].between(year_min, year_max) &
        #     d["Study"].isin(input.study_selector()) &
        #     d["Region"].isin(input.region_selector()) &
        #     (d["resistance_score"] >= min_resistance)
        # ]
        # 
        # if filtered_d.empty:
        #     ui.notification_show("No strains matching filter", type="warning")

        return d.data


#     @reactive.Effect
#     @reactive.event(reactive_loaded_data)
#     def _update_region_selector():
#         d = reactive_loaded_data()
#         if d is None or d.empty or "Region" not in d.columns:
#             regions = []
#         else:
#             regions = sorted(d["Region"].dropna().unique().tolist())
#
#         ui.update_select(session, "region_selector", choices=regions, selected=regions)
#
#     # Update year slider
#     @reactive.Effect
#     @reactive.event(reactive_loaded_data)
#     def _update_year_slider():
#         d = reactive_loaded_data()
#         if d is None or d.empty or "Year" not in d.columns:
#             min_y, max_y = 1900, 2100
#         else:
#             years = d["Year"].dropna()
#             min_y, max_y = int(years.min()), int(years.max())
#
#         ui.update_slider("year_selector", min=min_y, max=max_y, value=(min_y, max_y))
#
#     # Data selector reset ----------------------------------------------------------
#     @reactive.Effect
#     @reactive.event(input.data_reset)
#     def _reset_data_filters():
#         ui.notification_show("Resetting data filters", type="message")
#         ui.update_select(session, "amr_selector", selected="All samples")
#         ui.update_select(session, "year_variable", selected=[])
#         ui.update_select(session, "region_variable", selected=[])
#         ui.update_select(session, "country_variable", selected=[])
#
#         # Re-trigger updates for study, region, year
#         _update_study_selector()
#         _update_region_selector()
#         _update_year_slider()
#
#
#     # Control panels ---------------------------------------------------------------
    @reactive.effect
    @reactive.event(reactive_dataset)
    def _toggle_panels_on_load():
        d: Dataset | None = reactive_dataset.get()
        # TODO: Always remove panels first to handle re-uploads
        # ui.remove_ui(selector="prevalence_panel", immediate=True)
        # ui.remove_ui(selector="coverage_panel", immediate=True)
        # ui.remove_ui(selector="dataframe_panel", immediate=True)
        if d is None or len(d) == 0:
            ui.update_sidebar("sidebar", show=False)
            ui.update_accordion("upload_panel", show=False)
        else:
            ui.update_sidebar("sidebar", show=True)
            ui.update_accordion("upload_panel", show=True)
            metadata_columns = list(d.metadata_columns)
            genotype_columns = list(d.genotype_columns)
            ui.update_selectize('genotype_variable', choices=genotype_columns),
            ui.update_selectize('country_variable', choices=metadata_columns),
            ui.update_selectize('region_variable', choices=metadata_columns),

#
#     # Output dataframe -------------------------------------------------------------
#     @render.data_frame
#     def dataframe():
#         df = reactive_data()
#         if df is None or df.empty:
#             return
#         return render.DataGrid(df.reset_index(drop=True))
#
#     # Output summary ---------------------------------------------------------------
#     @render.text
#     def summary():
#         loaded = reactive_loaded_data()
#         if loaded is None or loaded.empty:
#             return ""
#
#         filtered = reactive_data()
#         if filtered is None or filtered.empty:
#             return f"Samples: 0/{len(loaded)}"
#
#         # This is a simplified version of your R summary string
#         return (
#             f"Samples: {len(filtered)}/{len(loaded)}; "
#             f"Studies: {len(input.study_selector() or [])}/{len(loaded['Study'].unique())}; "
#             f"Regions: {len(input.region_selector() or [])}/{len(loaded['Region'].unique())}; "
#             f"Years: {input.year_selector()[0]}-{input.year_selector()[1]}; "
#             f"Resistance: {input.amr_selector()}"
#         )
#
#     @reactive.calc
#     def reactive_data():
#         """Filters the main dataframe based on UI inputs."""
#         filtered_df = DATA[
#             (DATA['Year'].between(input.year_selector()[0], input.year_selector()[1])) &
#             (DATA['Study'].isin(input.study_selector())) &
#             (DATA['Region'].isin(input.region_selector())) &
#             (DATA['resistance_score'] >= RESMAP.get(input.amr_selector(), 0)) &
#             (DATA['Mortality'] == 1 if input.outcome_selector() != 'All samples' else True)
#             ]
#         if filtered_df.empty:
#             ui.notification_show("No strains matching filter", type="warning")
#         return filtered_df
#
#     @reactive.calc
#     def reactive_prevalence() -> PrevalenceResult:
#         """Calculates overall prevalence for the selected genotype."""
#         return PrevalenceCalculator(
#             stratify_by=[input.genotype_selector()],
#             adjust_for=['Cluster'],
#             n_distinct=[v for v in X_AXIS_VARIABLES if v not in {
#                 input.genotype_selector(), input.heatmap_x()}]
#         ).calculate(reactive_data())
#
#     @reactive.calc
#     def reactive_prevalence_stratified() -> PrevalenceResult:
#         """Calculates prevalence stratified by a second variable for the heatmap."""
#         return PrevalenceCalculator(
#             stratify_by=[input.genotype_selector(), input.heatmap_x()],
#             adjust_for=['Cluster'],
#             denominator=(
#                 input.genotype_selector()
#                 if input.heatmap_swap_denominator()
#                 else input.heatmap_x()
#             )
#         ).calculate(reactive_data())
#
#     @reactive.calc
#     def reactive_prevalence_coverage() -> PrevalenceResult:
#         """Calculates prevalence stratified by a second variable for the heatmap."""
#         return PrevalenceCalculator(
#             stratify_by=[input.coverage_group(), input.genotype_selector()],
#             adjust_for=['Cluster']
#         ).calculate(reactive_data())
#
#     @reactive.calc
#     def genotype_list() -> list[str]:
#         """Determines the list of genotypes to display based on user selections."""
#         if len(result := reactive_prevalence()) == 0:
#             return []
#         valency = min(input.vaccine_valency(), len(result))
#         return result.data[input.genotype_selector()].iloc[:valency].tolist()
#
#     @output
#     @render_widget
#     def global():
#         """Renders the main combined plot (pyramid, heatmap, bars)."""
#         # Get reactive variables
#         r1 = reactive_prevalence()
#         r2 = reactive_prevalence_stratified()
#         y_order = input.genotype_selected()
#         if len(r1) == 0 or len(r2) == 0 or not y_order:
#             return None
#         # Init plotters
#         p1 = PrevalencePlotter(y_order=y_order)
#         p2 = StrataPlotter(y_order=y_order, max_x=input.heatmap_num_x())
#         p3 = SummaryBarPlotter(y_order=y_order, fill_by=input.bars_x())
#         # Render plots
#         return merge_prevalence_figs(p1.plot(r1), p2.plot(r2), p3.plot(r1))
#
#     @output
#     @render_widget
#     def coverage_plot():
#         """Renders the main combined plot (pyramid, heatmap, bars)."""
#         if len(r := reactive_prevalence_coverage()) == 0 or not (x_order := genotype_list()):
#             return None
#         return CoveragePlotter(x_order=x_order).plot(r)
#
#     @output
#     @render_widget
#     def map():
#         """Renders the main combined plot (pyramid, heatmap, bars)."""
#         if len(r := reactive_prevalence_coverage()) == 0:
#             return None
#         return  MapPlotter().plot(r)
#
#     @output
#     @render.text
#     def summary():
#         """
#         Returns a summary of the filtered data
#         """
#         return (f"Samples: {len(reactive_data())}/{len(DATA)}; "
#             f"Studies: {len(input.study_selector())}/{len(DATA['Study'].unique())}; "
#             f"Regions: {len(input.region_selector())}/{len(DATA['Region'].unique())}; "
#             f"Years: {input.year_selector()[1] - input.year_selector()[0] + 1}/"
#             f"{DATA['Year'].max() - DATA['Year'].min() + 1}; "
#             f"Resistance: {input.amr_selector()}; "
#             f"Outcome: {input.outcome_selector()}")
