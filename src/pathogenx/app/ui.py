"""
Module for dealing with the web-app user interface
"""
from shiny import ui
from typing import get_args
from faicons import icon_svg as icon
from shinyswatch import theme
from pathogenx.app.utils import create_logo_link
from pathogenx.io import _GENOTYPE_FLAVOURS, _META_FLAVOURS, _DIST_FLAVOURS

_FLAVOURS = {'genotype': get_args(_GENOTYPE_FLAVOURS), 'metadata': get_args(_META_FLAVOURS), 'distance': get_args(_DIST_FLAVOURS)}
_VAR_CATEGORIES = ('genotype', 'adjustment', 'spatial', 'temporal', 'custom')


# Define all hyperlinked images here -------------------------------------------
kaptive_logo = create_logo_link("kaptive.png", "kaptive.readthedocs.io", "100px", 'Read the docs')
kleborate_logo = create_logo_link("kleborate.png", "kleborate.readthedocs.io", "100px", 'Read the docs')
monash_logo = create_logo_link("monash.svg", "www.monash.edu", "100px")
lshtm_logo = create_logo_link("lshtm.png", "www.lshtm.ac.uk", "70px")

# Define footer ----------------------------------------------------------------
footer = ui.div(
    ui.div(
        ui.div(
            lshtm_logo, monash_logo,
            style="display: flex; align-items: center; flex-wrap: nowrap;"
        ),
        ui.div(
            kaptive_logo, kleborate_logo,
            style="display: flex; align-items: center; flex-wrap: nowrap;"
        ),
        style="display: flex; justify-content: space-between; align-items: center;",
    ),
    style="margin: 10px auto; width: 100%;",
)

# Define home panel ------------------------------------------------------------
home = ui.nav_panel(
    'Home',
    ui.row(
        ui.column(
            10,
            ui.card(
                ui.card_body(
                    ui.HTML("""
                    <p>This app allows users to explore the distribution of predicted K  and O serotypes for 
                    <i>Klebsiella pneumoniae</i> isolated from neonatal sepsis cases in 13 studies across countries in 
                    Africa and Southern Asia, reported in the paper Stanton et al, 2025.</p>
                    <p>The functionality is geared towards exploring sets of K/O genotypes, in terms of their 
                    prevalence and distribution across geographical regions and theoretical coverage of infection 
                    isolates, to inform vaccine design.</p> 
                    <p>Prevalence estimates are adjusted for localised nosocomial clustering, to reduce the bias 
                    introduced by random outbreaks during surveillance periods. Coverage estimates are based on total 
                    isolates, not adjusted for clustering.</p>
                    <p>The <b>Modelled genotype prevalence</b> tab is populated with pre-calculated global and regional 
                    prevalence estimates modelled using Bayesian hierarchical meta-analysis, as described in the 
                    paper. Subgroup analyses are limited to those modelled and reported in the paper (geographic 
                    regions, fatal cases, ESBL- or carbapenemase- carrying isolates).</p>
                    <p>The <b>Dynamic genotype prevalence</b> tab is populated with simple pooled estimates calculated 
                    on the fly, allowing users to interactively explore prevalence and coverage more flexibly by 
                    country, study, year, and multi-locus sequence type (ST).</p>
                    """)
                )
            ),
            offset=1,
        ),
        style="margin-top: 20px;",
    ),
    ui.hr(),
    icon=icon('house')
)

# genotype selector -------------------------------------------------------------
def _setup_variables(i: str):
    return ui.input_selectize(f'{i}_variable', f'Select {i} variable', choices=[])

def _setup_filters(i: str):
    id_, label = f'{i}_filter', f'Filter by {i}'
    if i == 'temporal':
        return ui.input_slider(id_, label, min=1900, max=2050, value=(1900, 2050), step=1)
    return ui.input_selectize(id_, f'Filter by {i} variable', choices=[], multiple=True)

# Data filters -----------------------------------------------------------------
sidebar = ui.sidebar(
    ui.accordion(
        ui.accordion_panel(
            ui.h5("Variable selector"), *map(_setup_variables, _VAR_CATEGORIES),
            value="variable_selector", icon=icon("clock")
        ),
        ui.accordion_panel(
            ui.h5("Data filters"), *map(_setup_filters, _VAR_CATEGORIES),
            value="filter_selector", icon=icon("filter")
        ),
        id="sidebar_accordion",
        multiple=True,
        open="data_selector",
    ),
    width=350,
    id="sidebar",
    open="closed",
    # bg="white",
)

def _create_upload(ftype: str):
    return ui.input_file(f"{ftype}_upload", f"Select {ftype} file", accept=[".csv", ".txt", ".tsv"],
                         placeholder=f'{ftype}.csv')


def _create_flavours(i: tuple[str, list[str]]):
    ftype, flav = i
    return ui.input_select(f"{ftype}_flavour", f"{ftype[0].upper()}{ftype[1:]} file flavour", choices=flav,
                           selected=flav[0])

def _create_sub_panel(title, value=None, icon_=None, content: list=None):
    stem = title.split()[0].lower()
    value = value or f'{stem}_panel'
    content = content or [ui.output_ui(f"{value}_content")]
    return ui.accordion_panel(ui.h4(title), *content, value=value, icon=icon(icon_ or stem))

upload_panel = _create_sub_panel(
    "Upload files", content=[
        ui.row(
            ui.column(4, *map(_create_upload, _FLAVOURS.keys())),
            ui.column(4, *map(_create_flavours, _FLAVOURS.items())),
            ui.column(
                4,
                ui.input_slider("snp_distance", "SNP distance for clustering", min=0, max=100, value=20, step=1),
                ui.input_select("cluster_method", "Select clustering method", choices=['connected_components']),
                ui.input_action_button('load_data', 'Load data', class_='btn-primary', width='300px'),
                ui.hr(),
                ui.input_action_button('upload_reset', 'Reset uploads', class_='btn-danger', width='300px'),
            )
        ),
    ]
)

# Main panel -------------------------------------------------------------------
main_panel = ui.nav_panel(
    "Explore genotype prevalence and coverage",
    ui.layout_sidebar(
        sidebar,
        ui.HTML(
            "<p>Upload Kleborate, difference matrix, and metadata files from "
            "<a href='https://pathogen.watch/'>PathogenWatch</a> collections here! You can upload Kleborate CSV files "
            "from as many collections as you want, with optional difference matrix for cluster-adjusted prevalence, "
            "and metadata for spatio-temporal coverage breakdown. Files downloaded from "
            "PathogenWatch must be unmodified, i.e. remain in CSV format, and have "
            "unmodified names and content.</p>"
            "<p>Upon upload, valid files from each collection will be grouped together, "
            "available metadata will be added and clusters will be calculated from "
            "the difference matrix if provided. You will then have the option of filtering "
            "the data and selecting spatio-temporal variables for regional prevalence or "
            "further filtering.</p>"
        ),
        # ui.output_text("summary"),
        ui.accordion(
            upload_panel,
            _create_sub_panel("Total prevalence", "prevalence_panel", "earth-africa"),
            _create_sub_panel("Spatial coverage", "coverage_panel", "map"),
            _create_sub_panel("Table", "dataframe_panel", "table"),
            id="accordion", multiple=True, open=True
        ),
    ),
    icon=icon("laptop"),
)

# Define the main UI -----------------------------------------------------------
main_ui = ui.page_navbar(
    home,
    main_panel,
    title='PathoGenX',
    footer=footer,
    theme=theme.lumen,
    fillable=False,
    window_title="PathoGenX"
)
