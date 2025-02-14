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
                   + "(IllustrisTNG), and detector sensitivity (of O3, O4low, "\
                   + "O4high, design of three gravitational wave detectors "\
                   + "H1, L1, V1).", #TODO
    'md5': None, #TODO
    'title': "Auxiliary POSYDON data", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}

# v1
ZENODO_COLLECTION['v1_for_v2.0.0-pre1'] = {
    'data': "https://zenodo.org/record/14205146/files/POSYDON_data.tar.gz",
    'description': "The POSYDON v1 dataset post-processed to work with the "\
                   + "v2.0.0-pre1 version of the code. It contains "\
                   + "single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, CO-HeMS, "\
                   + "CO-HeMS_RLO all at solar metallicity. Additionally, it "\
                   + "includes data on supernova prescriptions "\
                   + "(Sukhbold+2016, Couch+2020, Patton+Sukhbold2020), star "\
                   + "formation history (IllustrisTNG), and detector "\
                   + "sensitivity (of O3, O4low, O4high, design of three "\
                   + "gravitational wave detectors H1, L1, V1).",
    'md5': "cf645a45b9b92c2ad01e759eb1950beb",
    'title': "Re-postprocessed POSYDON v1.0 dataset compatible with code "\
             + "release v2.0.0-pre1",
    'url': "https://zenodo.org/records/14205146"
}
ZENODO_COLLECTION['super-Eddington_v1'] = {
    'data': "https://zenodo.org/records/14216817/files/"\
            + "POSYDON_data_super_Eddington.tar.gz",
    'description': "Data for super-Eddington accretion compatible with "\
                   + "'v1_for_v2.0.0-pre1'. It contains CO-HMS_RLO, CO-HeMS, "\
                   + "CO-HeMS_RLO all at solar metallicity to replace the "\
                   + "corresponding data in 'v1_for_v2.0.0-pre1'",
    'md5': "64b46fbd65cb54819371a1cbc714019c",
    'title': "Re-postprocessed POSYDON v1.0 dataset, assuming "\
             + "super-Eddington accretion, compatible with code release "\
             + "v2.0.0-pre1",
    'url': "https://zenodo.org/records/14216817"
}
COMPLETE_SETS['v1'] = ['v1_for_v2.0.0-pre1']
COMPLETE_SETS['super-Eddington_v1'] = ['v1_for_v2.0.0-pre1',
                                       'super-Eddington_v1']

# v2
ZENODO_COLLECTION['v2_girds_2Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at twice solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at twice solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 2 Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_1Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_0.45Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 0.45 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 0.45 solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 0.45 Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_0.2Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 0.2 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 0.2 solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 0.2 Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_0.1Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 0.1 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 0.1 solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 0.1 Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_0.01Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 0.01 solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 0.01 solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 0.01 Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_1e-3Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 10^{-3} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 10^{-3} solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 10^{-3} Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
ZENODO_COLLECTION['v2_girds_1e-4Zsun'] = {
    'data': None, #TODO
    'description': "The POSYDON v2 dataset at 10^{-4} solar metallicity. It "\
                   + "contains single-HMS, single-HeMS, HMS-HMS, CO-HMS_RLO, "\
                   + "CO-HeMS, CO-HeMS_RLO all at 10^{-4} solar metallicity.", #TODO
    'md5': None, #TODO
    'title': "POSYDON v2.0 dataset at 10^{-4} Zsun", #TODO
    'url': "https://zenodo.org/communities/posydon" #TODO
}
COMPLETE_SETS['v2'] = ['auxiliary', 'v2_girds_2Zsun', 'v2_girds_1Zsun',
                       'v2_girds_0.45Zsun', 'v2_girds_0.2Zsun',
                       'v2_girds_0.1Zsun', 'v2_girds_0.01Zsun',
                       'v2_girds_1e-3Zsun', 'v2_girds_1e-4Zsun']
COMPLETE_SETS['v2_2Zsun'] = ['auxiliary', 'v2_girds_2Zsun']
COMPLETE_SETS['v2_1Zsun'] = ['auxiliary', 'v2_girds_1Zsun']
COMPLETE_SETS['v2_0.45Zsun'] = ['auxiliary', 'v2_girds_0.45Zsun']
COMPLETE_SETS['v2_0.2Zsun'] = ['auxiliary', 'v2_girds_0.2Zsun']
COMPLETE_SETS['v2_0.1Zsun'] = ['auxiliary', 'v2_girds_0.1Zsun']
COMPLETE_SETS['v2_0.01Zsun'] = ['auxiliary', 'v2_girds_0.01Zsun']
COMPLETE_SETS['v2_1e-3Zsun'] = ['auxiliary', 'v2_girds_1e-3Zsun']
COMPLETE_SETS['v2_1e-4Zsun'] = ['auxiliary', 'v2_girds_1e-4Zsun']
