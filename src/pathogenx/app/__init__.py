from pathlib import Path
from shiny import App
from pathogenx.app.ui import main_ui
from pathogenx.app.server import main_server

app = App(main_ui, main_server, static_assets=Path(__file__).parent / "www")

# from re import compile as regex
# import pandas as pd
#
# # Constants ------------------------------------------------------------------------------------------------------------
# # Define regexes ---------------------------------------------------------------
# LOCUS_TYPE_REGEX = regex(r'[KO][L]?[0-9/vOKabcfg⍺β𝛾,]{1,6}') # Misses -D1
# YEAR_REGEX = regex(r'\b(19|20)\d{2}\b') # Realistic year regex
#
# # Cosmetic globals -----------------------------------------------------------------
# DEFAULT_VALENCY = 20
# SIDEBAR_WIDTH = 350
# LO_FILL = '#2E91E5'
# HI_FILL = '#FD3216'
# BAR_FILL = '#2E91E5'
# SPINNER_COLOR = "#0dc5c1"
# SPINNER_TYPE = 7
# SPINNER_SIZE = 1
#
#
#
# # Data globals -----------------------------------------------------------------
# DATA = pd.read_table('src/klebnnsapp/data/TableS3_IsolatesIncluded_NEEDSACCESSIONS.tsv', sep='\t', index_col=0).drop(
#     ['neonatal', 'N50','contig_count', 'total_size', 'K_locus_confidence',  'O_locus_confidence',
#      'species',  'Cluster7', 'Cluster14', 'Cluster56', 'Cluster365', 'Cluster1SNP', 'Cluster2SNP', 'Cluster3SNP',
#      'Cluster5SNP', 'Cluster10SNP', 'Cluster15SNP', 'Cluster20SNP', 'Cluster25SNP', 'Cluster50SNP'], axis=1
# )
# DATA['Site'] = DATA['Site'].astype(str)
#
#
# # Variable globals -------------------------------------------------------------
# KLEBORATE_VARIABLES = ["ST", "K_locus", "K_type", "O_locus", "O_type", "resistance_score"]
# X_AXIS_VARIABLES = ['Country', 'Year', 'Region', 'Study'] + KLEBORATE_VARIABLES
# COUNTRIES = pd.unique(DATA['Country']).tolist()
# YEARS = pd.unique(DATA['Year']).tolist()
# MIN_YEARS = min(YEARS)
# MAX_YEARS = max(YEARS)
# REGIONS = pd.unique(DATA['Region']).tolist()
# STUDIES = pd.unique(DATA['Study']).tolist()
# SITES = pd.unique(DATA['Site']).tolist()
# ST = pd.unique(DATA['ST']).tolist()
# KL = pd.unique(DATA['K_locus']).tolist()
# OL = pd.unique(DATA['O_locus']).tolist()
# AMR = ['All samples', 'ESBL+', 'Carbapenemase+']
# RESMAP = {i: n for n, i in enumerate(AMR)}
# # MODELLED_INFECTION_TYPES <- unique(MODELLED$Infection)

