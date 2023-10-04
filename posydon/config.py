import os
from dotenv import load_dotenv

load_dotenv()

# POSYDON environment variables
PATH_TO_POSYDON = os.getenv("PATH_TO_POSYDON")
PATH_TO_POSYDON_DATA = os.getenv("PATH_TO_POSYDON_DATA")