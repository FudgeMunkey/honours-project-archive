import os
from dotenv import load_dotenv
from hypothesis import settings, Verbosity

load_dotenv()

all_profiles = os.getenv("ALL_PROFILES").split(" ")

for profile in all_profiles:
    settings.register_profile(
        profile, max_examples=int(os.getenv(f"{profile.upper()}_EXAMPLES"))
    )

settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "dev"))
