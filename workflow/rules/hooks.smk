# Adding handlers for displaying status of the 
# pipeline and for getting job information for 
# previously submitted jobs using `jobby`:
# https://github.com/OpenOmics/scribble/blob/main/scripts/jobby/jobby
if config['options']['mode'] == 'slurm':
    onstart:
        shell(
            """
            # Move any job information for a previous
            # instance of the pipeline to logfiles
            sleep 5; rm -f COMPLETED FAILED RUNNING;
            touch RUNNING
            for f in job_information_*.tsv; do
                # Skip over non-existant files
                [ -e "${{f}}" ] || continue
                mv ${{f}} logfiles/;
            done

            for f in failed_jobs_*.tsv; do
                # Skip over non-existant files
                [ -e "${{f}}" ] || continue
                mv ${{f}} logfiles/;
            done
            """
        )

    onsuccess:
        shell(
            """
            # Get job information on all
            # previously submitted jobs 
            sleep 15; rm -f COMPLETED FAILED RUNNING;
            timestamp=$(date +"%Y-%m-%d_%H-%M-%S");
            ./workflow/scripts/jobby \\
                $(grep --color=never "^Submitted .* external jobid" logfiles/snakemake.log \\
                    | awk '{{print $NF}}' \\
                    | sed "s/['.]//g" \\
                    | sort \\
                    | uniq \\
                    | tr "\\n" " "
                ) \\
            > job_information_${{timestamp}}.tsv

            # Get information on any child 
            # job(s) that may have failed 
            grep --color=never \\
                '^jobid\\|FAILED' \\
                job_information_${{timestamp}}.tsv \\
            > failed_jobs_${{timestamp}}.tsv
            touch COMPLETED  
            """
        )

    onerror:
        shell(
            """
            # Get job information on all
            # previously submitted jobs 
            sleep 15; rm -f COMPLETED FAILED RUNNING;
            timestamp=$(date +"%Y-%m-%d_%H-%M-%S");
            ./workflow/scripts/jobby \\
                $(grep --color=never "^Submitted .* external jobid" logfiles/snakemake.log \\
                    | awk '{{print $NF}}' \\
                    | sed "s/['.]//g" \\
                    | sort \\
                    | uniq \\
                    | tr "\\n" " "
                ) \\
            > job_information_${{timestamp}}.tsv

            # Get information on any child 
            # job(s) that may have failed 
            grep --color=never \\
                '^jobid\\|FAILED' \\
                job_information_${{timestamp}}.tsv \\
            > failed_jobs_${{timestamp}}.tsv
            touch FAILED
            """
        )
