"""
Module for dealing with the web-app server logic
"""
import pandas as pd
from scipy.sparse import coo_matrix
from pathlib import Path
from shiny import Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import render_widget
from pathogenx.io import GenotypeFile, MetaFile, DistFile
from pathogenx.dataset import Dataset
from pathogenx.calculators import PrevalenceCalculator, PrevalenceResult
from pathogenx.app.plotters import (PrevalencePlotter, StrataPlotter, SummaryBarPlotter, CoveragePlotter, MapPlotter,
                                 merge_prevalence_figs)

_VAR_CATEGORIES = ('genotype', 'adjustment', 'spatial', 'temporal', 'custom')

def main_server(input: Inputs, output: Outputs, session: Session):
    # Reactive container for user uploaded files to be loaded ----------------------
    reactive_dataset = reactive.Value[Dataset | None](None)

    def _load_genotypes() -> pd.DataFrame | None:
        file_infos = input.genotype_upload()
        if not file_infos:
            ui.notification_show("Genotype file is required.", type="error")
            return None
        f = file_infos[0]
        genotype_file = GenotypeFile.from_flavour(Path(f["datapath"]), input.genotype_flavour())
        try:
            data = genotype_file.load()
            ui.notification_show('Successfully loaded genotype file', type='message')
            return data
        except Exception as e:
            ui.notification_show(f"Error loading genotype file: {e}", type="error", duration=None)
            return None

    def _load_metadata() -> pd.DataFrame | None:
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
            ui.notification_show(f"Error loading metadata file: {e}", type="error", duration=None)
            return None

    def _load_distances() -> tuple[coo_matrix, list[str]] | None:
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
            ui.notification_show(f"Error loading distance file: {e}", type="error", duration=None)
            return None

    @reactive.effect
    @reactive.event(input.load_data, ignore_none=False)
    def _load_data_and_create_dataset():
        """
        This event runs when the user clicks the 'Load data' button.
        It loads all files, creates a Dataset object, calculates clusters,
        and sets the reactive value, triggering downstream updates.
        """
        genotypes = _load_genotypes()
        if genotypes is None or genotypes.empty:
            return
        dataset = Dataset(genotypes, _load_metadata(), _load_distances())
        if dataset.distances is not None:
            ui.notification_show('Calculating clusters...')
            dataset.calculate_clusters(method=input.cluster_method(), distance=input.snp_distance())
        reactive_dataset.set(dataset)

    # User data to be filtered and used for prevalences ----------------------------
    @reactive.calc
    def reactive_data() -> pd.DataFrame | None:
        if (d := reactive_dataset.get()) is None or len(d) == 0:
            return None

        filtered_data = d.data  # This creates a copy via the Dataset.data attribute method

        for var in _VAR_CATEGORIES:
            if (variable_col := input[f"{var}_variable"]()) and (filter_values := input[f"{var}_filter"]()):
                if var == 'temporal':
                    min_val, max_val = filter_values
                    filtered_data = filtered_data[filtered_data[variable_col].between(min_val, max_val)]
                else:
                    print(filter_values)
                    filtered_data = filtered_data[filtered_data[variable_col].isin(filter_values)]

        if filtered_data.empty:
            ui.notification_show("No data matches the current filter selection.", type="warning")
            return None

        return filtered_data

    @reactive.effect
    @reactive.event(reactive_dataset)
    def _toggle_panels_on_load():
        d: Dataset | None = reactive_dataset.get()
        if d is None or len(d) == 0:
            ui.update_sidebar("sidebar", show=False)
            ui.update_accordion("upload_panel", show=False)
        else:
            ui.update_sidebar("sidebar", show=True)
            ui.update_accordion("upload_panel", show=True)
            metadata_cols = list(d.metadata_columns) if d.metadata_columns is not None else []
            genotype_cols = list(d.genotype_columns)
            all_cols = sorted(genotype_cols + metadata_cols)
            adjust_cols = ['Cluster'] if d.distances is not None else []
            for var, cols in zip(_VAR_CATEGORIES, (genotype_cols, adjust_cols, metadata_cols, metadata_cols, all_cols)):
                # Add a blank choice to allow the input to be unselected
                ui.update_selectize(f"{var}_variable", choices=[''] + cols)

    def _create_event_lambda(var_name: str):
        """Function factory to correctly capture the loop variable for the lambda."""
        return lambda: input[f"{var_name}_variable"]()

    @reactive.effect
    # This effect is explicitly triggered when any of the variable selection dropdowns change.
    # We use a function factory (_create_event_lambda) to avoid the classic Python closure-in-a-loop issue.
    @reactive.event(*(_create_event_lambda(var) for var in _VAR_CATEGORIES), ignore_init=True)
    def _update_filter_selectors():
        """
        Populates filter controls based on the unique values in the columns
        the user has chosen in the variable selectors.
        """
        if (d := reactive_dataset.get()) is None or len(d) == 0:
            return

        for var in _VAR_CATEGORIES:

            # Only proceed if a column has been selected from the dropdown.
            # This check handles both None and empty string ""
            if selected_col := input[f"{var}_variable"]():
                col_data = d.data[selected_col]
                if var == 'temporal':
                    min_, max_ = int(col_data.min()), int(col_data.max())
                    ui.update_slider(f"{var}_filter", min=min_, max=max_, value=(min_, max_))
                else:
                    choices = sorted(col_data.dropna().unique().tolist())
                    ui.update_selectize(f"{var}_filter", choices=choices, selected=[])

    @output
    @render.ui
    def prevalence_panel_content():
        if reactive_dataset.get() is not None:
            return ui.output_plot("merged_plot")

    @output
    @render.ui
    def coverage_panel_content():
        if input.spatial_variable() and input.temporal_variable():
            return ui.layout_column_wrap(
                ui.card(ui.card_body(ui.output_plot("coverage_plot"), class_="p-0"), full_screen=True),
                ui.card(ui.card_body(ui.output_ui("map"), class_="p-0"), full_screen=True),
                width=1 / 2,
                height="400px",
            )

    @output
    @render.ui
    def dataframe_panel_content():
        if input.spatial_variable() and input.temporal_variable():
            return ui.output_data_frame("dataframe")


    # Output dataframe -------------------------------------------------------------
    @render.data_frame
    def dataframe():
        if (df := reactive_data()) is None or df.empty:
            return None
        return render.DataGrid(df.reset_index(drop=True))
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

    @reactive.calc
    def reactive_prevalence() -> PrevalenceResult | None:
        """Calculates overall prevalence for the selected genotype."""
        if (dataset := reactive_dataset()) is None or len(dataset) == 0:
            return None
        if (data := reactive_data()) is None or len(data) == 0:
            return None
        if (genotype := input.genotype_variable()) is None:
            return None
        adjust_by = input.adjustment_variable()
        return PrevalenceCalculator(
            stratify_by=[genotype],
            adjust_for=[adjust_by] if adjust_by else None,
            n_distinct=list(dataset.genotype_columns - {genotype, input.heatmap_x()})
        ).calculate(data)

    @reactive.calc
    def reactive_prevalence_stratified() -> PrevalenceResult | None:
        """Calculates prevalence stratified by a second variable for the heatmap."""
        if (data := reactive_data()) is None or len(data) == 0:
            return None
        if (genotype := input.genotype_variable()) is None:
            return None
        adjust_by = input.adjustment_variable()
        return PrevalenceCalculator(
            stratify_by=[input.genotype_variable(), input.heatmap_x()],
            adjust_for=[adjust_by] if adjust_by else None,
            denominator=(genotype if input.heatmap_swap_denominator() else input.heatmap_x())
        ).calculate(data)

    @reactive.calc
    def reactive_prevalence_coverage() -> PrevalenceResult | None:
        """Calculates prevalence stratified by a second variable for the heatmap."""
        if (data := reactive_data()) is None or len(data) == 0:
            return None
        if (genotype := input.genotype_variable()) is None:
            return None
        if (spatial := input.spatial_variable()) is None:
            return None
        adjust_by = input.adjustment_variable()
        return PrevalenceCalculator(
            stratify_by=[spatial, genotype],
            adjust_for=[adjust_by] if adjust_by else None
        ).calculate(data)

    @output
    @render_widget
    def merged_plot():
        """Renders the main combined plot (pyramid, heatmap, bars)."""
        # Get reactive variables
        if (r1 := reactive_prevalence()) is None or len(r1) == 0:
            return None
        if (r2 := reactive_prevalence_stratified()) is None or len(r2) == 0:
            return None
        # Init plotters
        p1 = PrevalencePlotter()
        p2 = StrataPlotter(max_x=input.heatmap_num_x())
        p3 = SummaryBarPlotter(fill_by=input.bars_x())
        # Render plots
        return merge_prevalence_figs(p1.plot(r1), p2.plot(r2), p3.plot(r1))

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
