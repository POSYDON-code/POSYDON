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
                          'url': "https://zenodo.org/communities/posydon"
                         }}
COMPLETE_SETS = {}

# auxiliary data
ZENODO_COLLECTION['auxiliary'] = {
    'data': None, #TODO
    'description': "Auxiliary data for POSYDON. It contains data on "\
                   + "supernova prescriptions (Sukhbold+2016, Couch+2020, "\
                   + "Patton+Sukhbold2020), star formation history "\
                   + "(IllustrisTNG, Zavala+2021, Chruslinska+2021), and "\
                   + "detector sensitivity (of O3, O4low, O4high, design of "\
                   + "three gravitational wave detectors H1, L1, V1). About "\
                   + "4GB on disk.",
    'md5': None, #TODO
    'title': "Auxiliary POSYDON data",
    'url': "https://zenodo.org/communities/posydon" #TODO
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
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at twice solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 2 Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_1Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_0.45Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 0.45 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 0.45 Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_0.2Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 0.2 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 0.2 Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_0.1Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 0.1 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 0.1 Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_0.01Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 0.01 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 0.01 Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_1e-3Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 10^{-3} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 10^{-3} Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['DR2_grids_1e-4Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON DR2 dataset at 10^{-4} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, HMS-HMS_RLO,"\
                   + " CO-HMS_RLO, CO-HeMS, CO-HeMS_RLO all at twice solar "\
                   + "metallicity.",
    'md5': None, #TODO
    'title': "POSYDON DR2 dataset at 10^{-4} Zsun",
    'url': "https://zenodo.org/communities/posydon" #TODO
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
