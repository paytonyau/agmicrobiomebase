# ENA upload process

`UpENA.sh` is a tool designed to simplify the process of uploading large FASTQ files to the ENA. This tool, created by [Micha Bayer](https://www.hutton.ac.uk/staff/micha-bayer), has significantly streamlined the data upload process, and for this, we express our gratitude for his valuable contribution to the scientific community.

## Usage

The script requires three arguments:

1. `username`: Your ENA account username.
2. `password txt file`: A file containing your ENA account password. Example .txt file: `PW.txt`
3. `txt file with paths to FASTQ files for transfer`: A file containing the full paths to the FASTQ files you want to upload. Each path should be on a new line. Example .txt.file: `input_config.txt`

Here is the syntax to run the script:


```./UpENA.sh <username> <password file> <file with paths to FASTQ files for transfer>```

For example, if your username is ‘user123’, your password is stored in ‘pass.txt’, and the paths to the FASTQ files are in ‘files.txt’, you would run:

```./UpENA.sh user123 pass.txt files.txt```


## Output
The script will output the start and end time of the file upload, as well as the status of each file transfer.


## Error Handling
If the incorrect number of arguments is provided, the script will output an error message and exit.