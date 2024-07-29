

# FASTQ files contain raw sequence data and quality scores. 

# Use the Bio.SeqIO module from Biopython to parse FASTQ files.

from Bio import SeqIO

# SeqIO.parse: Reads the FASTQ file.

def parse_fastq(file_path):
    """
    Parses a FASTQ file and prints the sequence ID, sequence, and quality scores for each record.

    Args:
        file_path (str): Path to the FASTQ file.
    """
    for record in SeqIO.parse(file_path, "fastq"):

        print("ID:", record.id)
        # record.id: The sequence identifier.

        print("Sequence:", record.seq)
        # record.seq: The DNA sequence.

        print("Quality Scores:", record.letter_annotations["phred_quality"])
        # record.letter_annotations["phred_quality"]: The quality scores.
        

# Example usage
# parse_fastq("sample.fastq")






# BAM/SAM files store aligned sequence data. 
# We can use the pysam library to parse BAM/SAM files.

import pysam

def parse_bam(file_path):
    """
    Parses a BAM/SAM file and prints the query name, sequence, and flag for each read.

    Args:
        file_path (str): Path to the BAM/SAM file.
    """
    bamfile = pysam.AlignmentFile(file_path, "rb")
    # pysam.AlignmentFile: Opens the BAM/SAM file.

    for read in bamfile:
        print("Query Name:", read.query_name)
        # read.query_name: The name of the read.

        print("Sequence:", read.query_sequence)
        # read.query_sequence: The aligned sequence.

        print("Flag:", read.flag)
        # read.flag: The alignment flag indicating the status of the read.

# Example usage
# parse_bam("sample.bam")






# VCF files store genetic variants. We can use the vcfpy library to parse VCF files.

import vcfpy

def parse_vcf(file_path):
    """
    Parses a VCF file and prints the chromosome, position, ID, reference, and alternate alleles for each record.

    Args:
        file_path (str): Path to the VCF file.
    """
    reader = vcfpy.Reader.from_path(file_path)
    # vcfpy.Reader.from_path: Reads the VCF file.

    for record in reader:
        print("Chromosome:", record.CHROM)
        # record.CHROM: The chromosome.

        print("Position:", record.POS)
        # record.POS: The position on the chromosome.

        print("ID:", record.ID)
        # record.ID: The variant ID.

        print("Reference Allele:", record.REF)
        # record.REF: The reference allele.

        print("Alternate Alleles:", record.ALT)
        # record.ALT: The alternate alleles.

# Example usage
# parse_vcf("sample.vcf")






#  GFF/GTF files contain genome annotations. We can use the gffutils library to parse these files.

import gffutils

def parse_gff(file_path):
    """
    Parses a GFF/GTF file and prints the feature ID, sequence ID, start, end, and strand for each feature.

    Args:
        file_path (str): Path to the GFF/GTF file.
    """
    db = gffutils.create_db(file_path, dbfn=':memory:')
    # gffutils.create_db: Creates an in-memory database from the GFF/GTF file.

    for feature in db.all_features():
        print("Feature ID:", feature.id)
        # feature.id: The feature ID.

        print("Sequence ID:", feature.seqid)
        # feature.seqid: The sequence ID (chromosome).

        print("Start:", feature.start)
        # feature.start: The start position of the feature.

        print("End:", feature.end)
        # feature.end: The end position of the feature.

        print("Strand:", feature.strand)
        # feature.strand: The strand (+ or -).

# Example usage
# parse_gff("sample.gff")




def main():
    while True:
        print("Select the type of file to parse:")
        print("1. FASTQ")
        print("2. BAM/SAM")
        print("3. VCF")
        print("4. GFF/GTF")
        print("5. Exit")
        
        choice = input("Enter your choice (1-5): ")
        
        if choice == '5':
            print("Exiting the program.")
            break
        
        file_path = input("Enter the path to the file: ")
        
        if choice == '1':
            parse_fastq(file_path)
        elif choice == '2':
            parse_bam(file_path)
        elif choice == '3':
            parse_vcf(file_path)
        elif choice == '4':
            parse_gff(file_path)
        else:
            print("Invalid choice. Please select a number between 1 and 5.")

if __name__ == "__main__":
    main()



from Bio import SeqIO
import pysam
import vcfpy
import gffutils
import os

def parse_fastq(file_path):
    try:
        for record in SeqIO.parse(file_path, "fastq"):
            print("ID:", record.id)
            print("Sequence:", record.seq)
            print("Quality Scores:", record.letter_annotations["phred_quality"])
    except Exception as e:
        print(f"Error parsing FASTQ file: {e}")

def parse_bam(file_path):
    try:
        bamfile = pysam.AlignmentFile(file_path, "rb")
        for read in bamfile:
            print("Query Name:", read.query_name)
            print("Sequence:", read.query_sequence)
            print("Flag:", read.flag)
    except Exception as e:
        print(f"Error parsing BAM/SAM file: {e}")

def parse_vcf(file_path):
    try:
        reader = vcfpy.Reader.from_path(file_path)
        for record in reader:
            print("Chromosome:", record.CHROM)
            print("Position:", record.POS)
            print("ID:", record.ID)
            print("Reference Allele:", record.REF)
            print("Alternate Alleles:", record.ALT)
    except Exception as e:
        print(f"Error parsing VCF file: {e}")

def parse_gff(file_path):
    try:
        db = gffutils.create_db(file_path, dbfn=':memory:')
        for feature in db.all_features():
            print("Feature ID:", feature.id)
            print("Sequence ID:", feature.seqid)
            print("Start:", feature.start)
            print("End:", feature.end)
            print("Strand:", feature.strand)
    except Exception as e:
        print(f"Error parsing GFF/GTF file: {e}")

def validate_file(file_path):
    if not os.path.isfile(file_path):
        print("File does not exist. Please check the path and try again.")
        return False
    return True

def main():
    while True:
        print("\nSelect the type of file to parse:")
        print("1. FASTQ")
        print("2. BAM/SAM")
        print("3. VCF")
        print("4. GFF/GTF")
        print("5. Exit")
        
        choice = input("Enter your choice (1-5): ")
        
        if choice == '5':
            print("Exiting the program.")
            break
        
        file_path = input("Enter the path to the file: ")
        
        if not validate_file(file_path):
            continue
        
        if choice == '1':
            parse_fastq(file_path)
        elif choice == '2':
            parse_bam(file_path)
        elif choice == '3':
            parse_vcf(file_path)
        elif choice == '4':
            parse_gff(file_path)
        else:
            print("Invalid choice. Please select a number between 1 and 5.")

if __name__ == "__main__":
    main()



# to enhance the program further,several improvements can be made, 
# including error handling, 
# file validation, 
# better user feedback, 
# logging, and 
# additional functionalities such as 
# filtering or data export


# Error Handling:
# - Added try-except blocks around the parsing functions to catch and report errors gracefully.

# File Validation:
# - Added a validate_file function to check if the file exists before attempting to parse it.

# User Feedback: 
# -Improved user feedback by confirming file existence and handling invalid choices.


# Modularity: 
# - Functions are neatly separated, making the code easier to maintain and extend.


# Logging:

#  - Add logging to record processing steps and errors.

# Filtering and Exporting:

# -  Implement options to filter data and export results to a file.

# GUI: 

# - Create a graphical user interface for easier use, possibly using a library like Tkinter or PyQt.

# Unit Tests:

#  - Write unit tests to ensure the functions work correctly with various input cases.

# We'll use Python's logging module to log processing steps and errors.

import logging
from Bio import SeqIO
import pysam
import vcfpy
import gffutils
import os

# Configure logging
logging.basicConfig(filename='genomic_data_parser.log', level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(message)s')

def parse_fastq(file_path):
    try:
        for record in SeqIO.parse(file_path, "fastq"):
            logging.info(f"Parsing record ID: {record.id}")
            print("ID:", record.id)
            print("Sequence:", record.seq)
            print("Quality Scores:", record.letter_annotations["phred_quality"])
    except Exception as e:
        logging.error(f"Error parsing FASTQ file: {e}")
        print(f"Error parsing FASTQ file: {e}")

def parse_bam(file_path):
    try:
        bamfile = pysam.AlignmentFile(file_path, "rb")
        for read in bamfile:
            logging.info(f"Parsing read query name: {read.query_name}")
            print("Query Name:", read.query_name)
            print("Sequence:", read.query_sequence)
            print("Flag:", read.flag)
    except Exception as e:
        logging.error(f"Error parsing BAM/SAM file: {e}")
        print(f"Error parsing BAM/SAM file: {e}")

def parse_vcf(file_path):
    try:
        reader = vcfpy.Reader.from_path(file_path)
        for record in reader:
            logging.info(f"Parsing VCF record ID: {record.ID}")
            print("Chromosome:", record.CHROM)
            print("Position:", record.POS)
            print("ID:", record.ID)
            print("Reference Allele:", record.REF)
            print("Alternate Alleles:", record.ALT)
    except Exception as e:
        logging.error(f"Error parsing VCF file: {e}")
        print(f"Error parsing VCF file: {e}")

def parse_gff(file_path):
    try:
        db = gffutils.create_db(file_path, dbfn=':memory:')
        for feature in db.all_features():
            logging.info(f"Parsing GFF feature ID: {feature.id}")
            print("Feature ID:", feature.id)
            print("Sequence ID:", feature.seqid)
            print("Start:", feature.start)
            print("End:", feature.end)
            print("Strand:", feature.strand)
    except Exception as e:
        logging.error(f"Error parsing GFF/GTF file: {e}")
        print(f"Error parsing GFF/GTF file: {e}")

def validate_file(file_path):
    if not os.path.isfile(file_path):
        logging.error("File does not exist")
        print("File does not exist. Please check the path and try again.")
        return False
    return True

def main():
    while True:
        print("\nSelect the type of file to parse:")
        print("1. FASTQ")
        print("2. BAM/SAM")
        print("3. VCF")
        print("4. GFF/GTF")
        print("5. Exit")
        
        choice = input("Enter your choice (1-5): ")
        
        if choice == '5':
            logging.info("Exiting the program")
            print("Exiting the program.")
            break
        
        file_path = input("Enter the path to the file: ")
        
        if not validate_file(file_path):
            continue
        
        if choice == '1':
            parse_fastq(file_path)
        elif choice == '2':
            parse_bam(file_path)
        elif choice == '3':
            parse_vcf(file_path)
        elif choice == '4':
            parse_gff(file_path)
        else:
            logging.warning("Invalid choice")
            print("Invalid choice. Please select a number between 1 and 5.")

if __name__ == "__main__":
    main()



# 2. Implementing Filtering and Exporting
# We'll add simple filtering for sequence length and the ability to export parsed data to a CSV file.


import csv

def parse_fastq(file_path, min_length=0):
    try:
        filtered_records = []
        for record in SeqIO.parse(file_path, "fastq"):
            if len(record.seq) >= min_length:
                logging.info(f"Parsing record ID: {record.id}")
                filtered_records.append({
                    "ID": record.id,
                    "Sequence": str(record.seq),
                    "Quality Scores": record.letter_annotations["phred_quality"]
                })
        return filtered_records
    except Exception as e:
        logging.error(f"Error parsing FASTQ file: {e}")
        print(f"Error parsing FASTQ file: {e}")

# Similar modifications will be done for parse_bam, parse_vcf, and parse_gff...

def export_to_csv(data, file_path):
    keys = data[0].keys()
    with open(file_path, 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames=keys)
        dict_writer.writeheader()
        dict_writer.writerows(data)
    logging.info(f"Data exported to {file_path}")

def main():
    while True:
        print("\nSelect the type of file to parse:")
        print("1. FASTQ")
        print("2. BAM/SAM")
        print("3. VCF")
        print("4. GFF/GTF")
        print("5. Exit")
        
        choice = input("Enter your choice (1-5): ")
        
        if choice == '5':
            logging.info("Exiting the program")
            print("Exiting the program.")
            break
        
        file_path = input("Enter the path to the file: ")
        
        if not validate_file(file_path):
            continue
        
        min_length = int(input("Enter minimum sequence length to filter (0 for no filter): "))
        export_path = input("Enter the path to export filtered data (leave blank to skip export): ")
        
        parsed_data = None
        if choice == '1':
            parsed_data = parse_fastq(file_path, min_length)
        elif choice == '2':
            parsed_data = parse_bam(file_path)
        elif choice == '3':
            parsed_data = parse_vcf(file_path)
        elif choice == '4':
            parsed_data = parse_gff(file_path)
        else:
            logging.warning("Invalid choice")
            print("Invalid choice. Please select a number between 1 and 5.")
        
        if parsed_data and export_path:
            export_to_csv(parsed_data, export_path)

if __name__ == "__main__":
    main()


# 3. Creating a GUI
# We'll use Tkinter to create a simple graphical user interface.


import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import simpledialog

def browse_file():
    file_path = filedialog.askopenfilename()
    return file_path

def process_file(file_type, min_length, export_path):
    file_path = browse_file()
    if not validate_file(file_path):
        messagebox.showerror("Error", "File does not exist.")
        return

    parsed_data = None
    if file_type == "FASTQ":
        parsed_data = parse_fastq(file_path, min_length)
    elif file_type == "BAM/SAM":
        parsed_data = parse_bam(file_path)
    elif file_type == "VCF":
        parsed_data = parse_vcf(file_path)
    elif file_type == "GFF/GTF":
        parsed_data = parse_gff(file_path)
    
    if parsed_data and export_path:
        export_to_csv(parsed_data, export_path)
        messagebox.showinfo("Success", f"Data exported to {export_path}")

def main_gui():
    root = tk.Tk()
    root.title("Genomic Data Parser")

    tk.Label(root, text="Select file type:").pack()
    file_type_var = tk.StringVar(value="FASTQ")
    tk.Radiobutton(root, text="FASTQ", variable=file_type_var, value="FASTQ").pack(anchor=tk.W)
    tk.Radiobutton(root, text="BAM/SAM", variable=file_type_var, value="BAM/SAM").pack(anchor=tk.W)
    tk.Radiobutton(root, text="VCF", variable=file_type_var, value="VCF").pack(anchor=tk.W)
    tk.Radiobutton(root, text="GFF/GTF", variable=file_type_var, value="GFF/GTF").pack(anchor=tk.W)
    
    tk.Label(root, text="Minimum sequence length (0 for no filter):").pack()
    min_length_var = tk.IntVar(value=0)
    tk.Entry(root, textvariable=min_length_var).pack()
    
    tk.Label(root, text="Export path (leave blank to skip export):").pack()
    export_path_var = tk.StringVar()
    tk.Entry(root, textvariable=export_path_var).pack()
    
    tk.Button(root, text="Process", command=lambda: process_file(
        file_type_var.get(), min_length_var.get(), export_path_var.get())).pack()

    root.mainloop()

if __name__ == "__main__":
    main_gui()


# 4. will use the unittest framework to write unit tests for our parsing functions.

# We'll write unit tests for each of the parsing functions to ensure they handle various input cases correctly.

# First, we need some sample data files for testing. We'll assume that these sample files are stored in a directory named tests.

# Sample Test Data Files
# tests/sample.fastq
# tests/sample.bam
# tests/sample.vcf
# tests/sample.gff
# Now, let's write the unit tests.

# test_genomic_data_parser.py

import unittest
from genomic_data_parser import parse_fastq, parse_bam, parse_vcf, parse_gff

class TestGenomicDataParser(unittest.TestCase):

    def test_parse_fastq(self):
        records = parse_fastq('tests/sample.fastq')
        self.assertIsInstance(records, list)
        self.assertGreater(len(records), 0)
        self.assertIn('ID', records[0])
        self.assertIn('Sequence', records[0])
        self.assertIn('Quality Scores', records[0])

    def test_parse_bam(self):
        records = parse_bam('tests/sample.bam')
        self.assertIsInstance(records, list)
        self.assertGreater(len(records), 0)
        self.assertIn('Query Name', records[0])
        self.assertIn('Sequence', records[0])
        self.assertIn('Flag', records[0])

    def test_parse_vcf(self):
        records = parse_vcf('tests/sample.vcf')
        self.assertIsInstance(records, list)
        self.assertGreater(len(records), 0)
        self.assertIn('Chromosome', records[0])
        self.assertIn('Position', records[0])
        self.assertIn('ID', records[0])
        self.assertIn('Reference Allele', records[0])
        self.assertIn('Alternate Alleles', records[0])

    def test_parse_gff(self):
        records = parse_gff('tests/sample.gff')
        self.assertIsInstance(records, list)
        self.assertGreater(len(records), 0)
        self.assertIn('Feature ID', records[0])
        self.assertIn('Sequence ID', records[0])
        self.assertIn('Start', records[0])
        self.assertIn('End', records[0])
        self.assertIn('Strand', records[0])

if __name__ == '__main__':
    unittest.main()
