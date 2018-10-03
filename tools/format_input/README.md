## Using these model format conversion tools
These Python tools take a model input and re-format it as a sedona input.
For a hypothetical code called "Great Code", the script to turn Great Code
output into Sedona input is `great_code.py`.

These conversion scripts use command line arguments to manage the variety
of output files possible from the origin code. 
To see the options possible, type, e.g.
`python great_code.py -h`
To run the conversion script:
`python great_code.py <flags> <filename>`

## Developing a new conversion script
These scripts should be compatible with Python 3. 

