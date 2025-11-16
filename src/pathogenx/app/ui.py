from shiny import ui
from faicons import icon_svg as icon

from pathogenx.app.utils import create_logo_link


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
                         <p>This app allows users to explore the distribution of predicted K 
                         and O serotypes for <i>Klebsiella pneumoniae</i> isolated from neonatal 
                         sepsis cases in 13 studies across countries in Africa and Southern 
                         Asia, reported in the paper Stanton et al, 2025.</p>

                         <p>The functionality is geared towards exploring sets of K/O genotypes, 
                         in terms of their prevalence and distribution across geographical 
                         regions and theoretical coverage of infection isolates, to inform 
                         vaccine design.</p> 

                         <p>Prevalence estimates are adjusted for localised nosocomial 
                         clustering, to reduce the bias introduced by random outbreaks during 
                         surveillance periods. Coverage estimates are based on total isolates, 
                         not adjusted for clustering.</p>

                         <p>The <b>Modelled genotype prevalence</b> tab is populated with 
                         pre-calculated global and regional prevalence estimates modelled 
                         using Bayesian hierarchical meta-analysis, as described in the 
                         paper. Subgroup analyses are limited to those modelled and reported 
                         in the paper (geographic regions, fatal cases, ESBL- or carbapenemase- carrying isolates).
                         </p>

                         <p>The <b>Dynamic genotype prevalence</b> tab is populated with simple pooled 
                         estimates calculated on the fly, allowing users to interactively explore prevalence 
                         and coverage more flexibly by country, study, year, and multi-locus sequence type 
                         (ST).</p>
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
genotype_selector = ui.accordion_panel(
    ui.h5("genotype selector"),
    ui.row(
        ui.column(
            6,
            ui.input_select(
                "genotype_selector",
                ui.h6("Select genotype:"),
                choices=["K_locus", "O_type"],
                selected="K_locus",
            ),
        ),
        ui.column(
            6,
            ui.input_slider(
                "vaccine_valency",
                ui.h6("# of genotypes"),
                min=1,
                max=150,  # Will be updated dynamically
                value=20,
            ),
        ),
    ),
    ui.p("Hit update to order genotypes by the most prevalent in the current data:"),
    ui.row(
        ui.column(
            10,
            ui.input_action_button(
                "genotype_reset",
                "Update genotype order",
                class_="btn-primary",
                icon=icon("spinner"),
                width="100%",
            ),
        ),
        ui.column(1, ui.output_ui("genotype_copy")),
    ),
    value="genotype_selector",
    icon=icon("syringe"),
)

# Data filters -----------------------------------------------------------------
sidebar = ui.sidebar(
    ui.accordion(
        genotype_selector,
        data_selector,
        id="sidebar_accordion",
        multiple=False,
        open="data_selector",
    ),
    width=350,
    id="sidebar",
    open="closed",
    bg="white",
)

upload_panel = ui.accordion_panel(
    ui.h4("Upload files"),
    ui.row(
        ui.column(
            8,
            ui.input_file(
                "genotype_upload", "Select genotype file",
                accept=[".csv", ".txt", ".tsv"], placeholder='kleborate.csv'
            ),
        ),
        ui.column(
            4,
            ui.input_select(
                "genotype_flavour", "Genotype file type",
                choices=['pw-kleborate', 'kleborate'], selected='pw-kleborate'
            )
        ),
    ),
    ui.row(
        ui.column(
            8,
            ui.input_file(
                "metadata_upload", "Select metadata file",
                accept=[".csv", ".txt", ".tsv"], placeholder='metadata.csv'
            ),
        ),
        ui.column(
            4,
            ui.input_select(
                "metadata_flavour", "Metadata file type",
                choices=['pw-metadata'], selected='pw-metadata'
            )
        ),
    ),
    ui.row(
        ui.column(
            8,
            ui.input_file(
                "distance_upload", "Select distance file",
                accept=[".csv", ".txt", ".tsv"], placeholder='distance.csv'
            ),
        ),
        ui.column(
            4,
            ui.input_select(
                "distance_flavour", "Genotype file type",
                choices=['pw-kleborate', 'kleborate'], selected='pw-dist'
            )
        ),
    ),
    value="upload_panel",
    icon=icon("upload"),
)

# ui.input_slider(
#                 "snvs", "SNV distance for clustering", min=0, max=100, value=10, step=5
#             ),

# Panels to be inserted dynamically
prevalence_panel = ui.accordion_panel(
    ui.h4("Total prevalence"),
    ui.output_plot("merged_plot", height="600px"),
    value="prevalence_panel",
    icon=icon("earth-africa"),
)

coverage_panel = ui.accordion_panel(
    ui.h4("Geographical coverage"),
    ui.input_select("coverage_group", "Calculate coverage of", choices=[], selected=None),
    ui.layout_column_wrap(
        1 / 2,
        ui.card(
            ui.card_body(
                ui.output_plot("coverage_plot"), class_="p-0"
            ),
            full_screen=True,
        ),
        ui.card(
            ui.card_body(
                ui.output_ui("map"), class_="p-0"
            ),
            full_screen=True,
        ),
        height="400px",
    ),
    value="coverage_panel",
    icon=icon("map")
)

dataframe_panel = ui.accordion_panel(
    ui.h4("Table"),
    ui.output_data_frame("dataframe"),
    value="dataframe_panel",
    icon=icon("table"),
)


# Main panel -------------------------------------------------------------------
main_panel = ui.nav_panel(
    "Analyse your own data",
    ui.layout_sidebar(
        sidebar,
        ui.p(
            "Upload Kleborate, difference matrix, and metadata files from ",
            ui.a("PathogenWatch", href="https://pathogen.watch/"),
            " collections here! You can upload Kleborate CSV files from as many collections "
            "as you want, with optional difference matrix for cluster-adjusted prevalence, "
            "and metadata for spatio-temporal coverage breakdown. Files downloaded from "
            "PathogenWatch must be unmodified, i.e. remain in CSV format, and have "
            "unmodified names and content.",
        ),
        ui.p(
            "Upon upload, valid files from each collection will be grouped together, "
            "available metadata will be added and clusters will be calculated from "
            "the difference matrix if provided. You will then have the option of filtering "
            "the data and selecting spatio-temporal variables for regional prevalences or "
            "further filtering."
        ),
        ui.output_text("summary"),
        ui.accordion(upload_panel, id="accordion", multiple=False, open=True),
    ),
    icon=icon("laptop"),
)

# Define the main UI -----------------------------------------------------------
main_ui = ui.page_navbar(
    home,
    main_panel,
    title='KlebSeroEpi',
    footer=footer,
    theme=ui.Theme(preset='pulse'),
    fillable=False,
    window_title="KlebSeroEpi"
)