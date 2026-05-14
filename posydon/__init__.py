from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("posydon")
except PackageNotFoundError:
    # Package is not installed
    __version__ = "unknown"

__author__ = "Tassos Fragos <Anastasios.Fragkos@unige.ch>"
__credits__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Nam Tran <tranhn03@gmail.com>",
    "Jaime Roman Garza <Jaime.Roman@etu.unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Juanga Serra Perez <jgserra@northwestern.edu>",
    "Philipp Moura Srivastava <philipp.msrivastava@gmail.com>",
    "Ying Qin <<yingqin2013@hotmail.com>",
    "Aaron Dotter <aaron.dotter@gmail.com>",
]
