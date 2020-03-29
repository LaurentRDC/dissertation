# Wrapper script to ensure that the conda 
# environment "dissertation" is always activated
activate dissertation | Out-Null
python .\build.py $args