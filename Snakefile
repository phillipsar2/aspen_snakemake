configfile: "path/to/config.yaml"

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        ""


# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"
