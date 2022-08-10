#!/usr/bin/env python3

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from sigProfilerPlotting import sample_portrait as sP
import matplotlib.pyplot as plt
import argparse, time, os


def sigProfilerPlot(input_folder, output_path, output_name, reference):
    print('Running SigProfiler with...', input_folder, output_path, output_name, reference)
    print('Creating SigProfiler Matrices')
    # Generate SigProfilerMatrix
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        output_name, 
        reference, 
        input_folder, 
        plot=False
    )
    print('Creating SigProfiler Sample Plots')
    # Create sample portrait plot
    try:
        # Fails with subsampled dataset 
        # or with samples that have low 
        # mutational burdens  
        sP.samplePortrait(input_folder, output_path, output_name, percentage=False)
    except Exception as e:
        print('WARNING: samplePortrait plotting function failed.',)
        ofh = "{}sample_portrait_{}.pdf".format(output_path, output_name)
        # Creating an empty pdf with warning message
        # Give time for previous file device to close
        time.sleep(5)
        f = plt.figure()
        plt.suptitle('WARNING: SigProfiler failed to create a sample portrait plot!')
        plt.title('{}'.format(os.path.basename(ofh)))
        plt.axis('off')
        f.savefig(ofh, bbox_inches='tight')


def checkRef(ref):
    ref = ref.lower()
    if ref == "hg19" or ref == "grch37":
        return "GRCh37"
    elif ref == "hg38" or ref == "grch38":
        return "GRCh38"
    elif ref == "mm10":
        return "mm10"
    else:
        raise argparse.ArgumentTypeError("Invalid reference genome")


def main(raw_args = None):
    # Parse args
    parser = argparse.ArgumentParser(description = "Run SigProfiler")

    # Reference genome
    parser.add_argument(
        "-r", "--ref",
        required=False,
        default='GRCh38',
        metavar="reference", 
        dest="ref",
        action="store", 
        type=str, 
        help="Reference Genome (GRCh37, GRCh38, mm10)"
    )

    # Directory with input files
    parser.add_argument(
        "-i", "--input", 
        metavar="input_folder",
        required=True,
        dest="input",
        action="store", 
        type=str, 
        help="Input folder containing vcf files"
    )

    # Output directory
    parser.add_argument(
        "-o", "--output", 
        metavar="output_path",
        required=True,
        dest="output_path",
        action="store", 
        type=str, 
        help="Path to save output"
    )

    # Used for output file creation
    parser.add_argument(
        "-p", "--project", 
        metavar="output_name",
        required=True,
        dest="output_name",
        action="store", 
        type=str, 
        help="Name to use when saving output file (project name)"
    )

    args = parser.parse_args(raw_args)

    # Create sigProfiler Plot
    sigProfilerPlot(
        args.input, 
        args.output_path, 
        args.output_name, 
        checkRef(args.ref)
    )


if __name__=='__main__':
    main()

