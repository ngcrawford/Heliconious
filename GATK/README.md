# Pipeline: #


    1.) Add read groups to each of the filtered/merged bam files. 
    2.) Merge bam files as necessary
    3.) Run GATK on subset first
        - I think
    
    More here...


# Scripts: #

*AddRGToBams.py* 

    usage: addRGToBAMs.py [-h] -i INPUT_DIR -l SAMPLE_LIBRAY

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_DIR, --input-dir INPUT_DIR
                            The input directory containing the bam files.
      -l SAMPLE_LIBRAY, --sample-libray SAMPLE_LIBRAY
                            The file containg the sample to run info.


