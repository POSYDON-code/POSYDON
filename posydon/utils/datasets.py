"""Collection of all POSYDON datasets on Zenodo

Each publication on Zenodo should have and entry in 'ZENODO_COLLECTION' with
    - 'data' : str or None
        The url to the data (archive) to download.
    - 'description' : str
        A informative help text to tell about the dataset.
    - 'md5' : str or None
        The md5 checksum corresponding to the file given in 'data'.
    - 'title' : str
        The title of the dataset.
    - 'url' : str
        The url to the Zenodo page of the dataset.

Individual datasets can be combined to a complete set as input for the POSYDON
code. The complete sets are collected as lists of individual sets, which get
layered in definition order.

"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

ZENODO_COLLECTION = {'POSYDON':
                         {'data': None,
                          'description': "The collection of all POSYDON data "\
                                         + "sets.",
                          'md5': None,
                          'title': "POSYDON community",
                          'url': "https://zenodo.org/records/15194708"
                         }}
COMPLETE_SETS = {}

# auxiliary data
ZENODO_COLLECTION['auxiliary'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_auxiliary.tar.gz",
    'description': "Auxiliary data for POSYDON. It contains data on "\
                   + "supernova prescriptions (Sukhbold+2016, Couch+2020, "\
                   + "Patton+Sukhbold2020), star formation history "\
                   + "(IllustrisTNG, Zavala+2021, Chruslinska+2021), and "\
                   + "detector sensitivity (of O3, O4low, O4high, design of "\
                   + "three gravitational wave detectors H1, L1, V1). About "\
                   + "4GB on disk.",
    'md5': "7eed35f08065628656fd140c03800cd6",
    'title': "Auxiliary POSYDON data",
    'url': "https://zenodo.org/records/15194708",
}

# DR1
ZENODO_COLLECTION['DR1_for_v2.0.0-pre1'] = {
    'data': "https://zenodo.org/record/14205146/files/POSYDON_data.tar.gz",
    'description': "The POSYDON DR1 dataset post-processed to work with the "\
                   + "v2.0.0-pre1 version of the code. It contains "\
                   + "single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, CO-HeMS, "\
                   + "CO-HeMS_RLO all at solar metallicity. Additionally, it "\
                   + "includes data on supernova prescriptions "\
                   + "(Sukhbold+2016, Couch+2020, Patton+Sukhbold2020), star "\
                   + "formation history (IllustrisTNG), and detector "\
                   + "sensitivity (of O3, O4low, O4high, design of three "\
                   + "gravitational wave detectors H1, L1, V1).",
    'md5': "cf645a45b9b92c2ad01e759eb1950beb",
    'title': "Re-postprocessed POSYDON DR1 dataset compatible with code "\
             + "release v2.0.0-pre1 and v2.0.0-pre2", #TODO create a version for code version v2.0.0 and update the entry
    'url': "https://zenodo.org/records/14205146"
}
ZENODO_COLLECTION['DR1-super_Eddington'] = {
    'data': "https://zenodo.org/records/14216817/files/"\
            + "POSYDON_data_super_Eddington.tar.gz",
    'description': "Data for super-Eddington accretion compatible with "\
                   + "'DR1_for_v2.0.0-pre1'. It contains CO-HMS_RLO, CO-HeMS, "\
                   + "CO-HeMS_RLO all at solar metallicity to replace the "\
                   + "corresponding data in 'DR1_for_v2.0.0-pre1'",
    'md5': "64b46fbd65cb54819371a1cbc714019c",
    'title': "Re-postprocessed POSYDON DR1 dataset, assuming "\
             + "super-Eddington accretion, compatible with code release "\
             + "v2.0.0-pre1 and v2.0.0-pre2", #TODO create a version for code version v2.0.0 and update the entry
    'url': "https://zenodo.org/records/14216817"
}
COMPLETE_SETS['DR1'] = ['DR1_for_v2.0.0-pre1']
COMPLETE_SETS['DR1.1'] = ['DR1_for_v2.0.0-pre1', 'DR1-super_Eddington']

# DR2
ZENODO_COLLECTION['DR2_grids_2Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_2Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at twice solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': "32d216cd95857b3d29c73188377f31dc",
    'title': "POSYDON DR2 dataset at 2 Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_1Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_1Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at solar "\
                   + "metallicity.",
    'md5': "fc6a17c9675e44909e3fae14aad4c4dc",
    'title': "POSYDON DR2 dataset at Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_0.45Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_0.45Zsun.tar.gz",
    'description': "The POSYDON v2 dataset at 0.45 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 0.45 solar metallicity.",
    'md5': "b2dc618a8399084c9c239b6c9af67bd7",
    'title': "POSYDON DR2 dataset at 0.45 Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_0.2Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_0.2Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at 0.2 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 0.2 solar "\
                   + "metallicity.",
    'md5': "37b6879d194f5ac5c0eb76b46330cfd8",
    'title': "POSYDON DR2 dataset at 0.2 Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_0.1Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_0.1Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at 0.1 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 0.1 solar "\
                   + "metallicity.",
    'md5': "07225092f2c21519c884afa33d143b59",
    'title': "POSYDON DR2 dataset at 0.1 Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_0.01Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_0.01Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at 0.01 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 0.01 solar "\
                   + "metallicity.",
    'md5': "25157c0e156c729b3a4cecbe878cbf51",
    'title': "POSYDON DR2 dataset at 0.01 Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_1e-3Zsun'] = {
    'data':"https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_1e-3Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at 10^{-3} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 10^{-3} solar "\
                   + "metallicity.",
    'md5': "6b5c260ac80a6abee20f9f5313fb3b83",
    'title': "POSYDON DR2 dataset at 10^{-3} Zsun",
    'url': "https://zenodo.org/records/15194708",
}
ZENODO_COLLECTION['DR2_grids_1e-4Zsun'] = {
    'data': "https://zenodo.org/records/15194708/files/POSYDON_data_v2_grids_1e-4Zsun.tar.gz",
    'description': "The POSYDON DR2 dataset at 10^{-4} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at 10^{-4} solar "\
                   + "metallicity.",
    'md5': "f44d0035d80ccc30998537de99d440e4",
    'title': "POSYDON DR2 dataset at 10^{-4} Zsun",
    'url': "https://zenodo.org/records/15194708",
}
COMPLETE_SETS['DR2'] = ['auxiliary', 'DR2_grids_2Zsun', 'DR2_grids_1Zsun',
                        'DR2_grids_0.45Zsun', 'DR2_grids_0.2Zsun',
                        'DR2_grids_0.1Zsun', 'DR2_grids_0.01Zsun',
                        'DR2_grids_1e-3Zsun', 'DR2_grids_1e-4Zsun']
COMPLETE_SETS['DR2_2Zsun'] = ['auxiliary', 'DR2_grids_2Zsun']
COMPLETE_SETS['DR2_1Zsun'] = ['auxiliary', 'DR2_grids_1Zsun']
COMPLETE_SETS['DR2_0.45Zsun'] = ['auxiliary', 'DR2_grids_0.45Zsun']
COMPLETE_SETS['DR2_0.2Zsun'] = ['auxiliary', 'DR2_grids_0.2Zsun']
COMPLETE_SETS['DR2_0.1Zsun'] = ['auxiliary', 'DR2_grids_0.1Zsun']
COMPLETE_SETS['DR2_0.01Zsun'] = ['auxiliary', 'DR2_grids_0.01Zsun']
COMPLETE_SETS['DR2_1e-3Zsun'] = ['auxiliary', 'DR2_grids_1e-3Zsun']
COMPLETE_SETS['DR2_1e-4Zsun'] = ['auxiliary', 'DR2_grids_1e-4Zsun']

ZENODO_COLLECTION['v2_tutorial_populations'] = {
    'data': "https://zenodo.org/records/15708476/files/POSYDON_tutorial_populations.tar.gz",
    'description': "The POSYDON v2 tutorial population data. It contains "\
                   + "example binary populations that can be used in the "\
                   + "tutorials of the code. It contians a 10.000 binary "\
                   + "population at each of the DR2 grid metallicities.",
    'md5': "d66d181bbd74d3cd30b4360729034566",
    'title': "POSYDON v2.0 tutorial population data",
    'url': "https://zenodo.org/records/15708476",
}
COMPLETE_SETS['v2_tutorial_populations'] = ['v2_tutorial_populations']

ZENODO_COLLECTION['2025_school_data'] = {
    'data': "https://zenodo.org/records/17902460/files/POSYDON_2025_school_data.tar.gz",
    'description': "The POSYDON School 2025 lab datasets. "\
                   + "These are various populations required to run the labs "\
                   + "from the school. Populations ranging from 10 to "\
                   + "1 million binaries are included, along with some small "\
                   + "MESA grid data for one of the labs.",
    'md5': "b5ee19b546a3377101efe61738902b70",
    'title': "POSYDON School 2025 population and MESA data",
    'url': "https://zenodo.org/records/17902460"
}
COMPLETE_SETS['2025_school_data'] = ['2025_school_data']