#!/bin/csh

set maxNJobs=30
set pause=120

set dir=$1
set out=$2
set api=$3
rm -rf ${out}/jpred_tmp
mkdir ${out}/jpred_tmp
# rm -rf ${dir}_error
# mkdir ${dir}_error
# mkdir ${dir}_output

set zero = 0
set counter=0
set nJobsRunning = 0
foreach line (`ls $dir/*.fasta`)
  @ counter = $counter + 1
  # sleep 2
  set lengthCheck = `cat $line | sed '1d' | wc -m`
  if ($lengthCheck<20 || $lengthCheck>800) then
    rm $line
    echo "  Length of sequence is outside of the allowed interval [20,800]: " $line
  else
    while ($nJobsRunning >= $maxNJobs)
	    #---------------------------------------------------------
	    sleep $pause
	    foreach lineid (`ls $dir/*.job_id`)
    		set id = `cat $lineid | awk '{print $6}'`
    		set name = `echo $lineid | sed 's%/% %' | awk '{print $2}' | sed 's%.fasta.job_id%%g' | awk '{print $1}'`
    		echo "  Going to check status of job: " $name
    		perl $api status jobid=$id getResults=no checkEvery=once silent > $lineid.log_status
    		set mystatus2 = `grep "ERROR" $lineid.log_status | wc -m`
    		if ($mystatus2 > $zero) then
    		    echo "Jod $id error-ed."
            rm -f $dir/$name.fasta $lineid*
    		    # mv $dir/$name.fasta ${dir}_error/.
    		    # mv -f $lineid* ${dir}_error/.
    		else
		      set mystatus = `grep "finished. Results available at the following" $lineid.log_status | wc -m`
          if ($mystatus > $zero) then
        		echo "  Jod $name finished. Getting results, creating job directory and moving there all the relevant files."
        		set myurl = `grep "www.compbio.dundee.ac.uk" $lineid.log_status`
            set myhtml = `echo $myurl | sed 's%results.html%concise%'`
            curl -s $myhtml > ${out}/jpred_tmp/$name
            rm -f $dir/$name.fasta $lineid*
            # or mv $lineid* ${dir}_output/.
          endif
        endif
	    end
	    #---------------------------------------------------------
	    set nJobsRunning = `ls $dir/*job_id | wc -l`
	  end
	  echo "  Input file: " $line "submitting..."
	  perl $api submit mode=single format=fasta silent longtime=on file=$line > $line.job_id
  endif
  set nJobsRunning = `ls $dir/*job_id | wc -l`
end

set moreLogs = 1
while ($moreLogs > 0)
	#---------------------------------------------------------
	sleep $pause
	foreach lineid (`ls $dir/*.job_id`)
    set id = `cat $lineid | awk '{print $6}'`
    set name = `echo $lineid | sed 's%/% %' | awk '{print $2}' | sed 's%.fasta.job_id%%g' | awk '{print $1}'`
    echo "  Going to check status of job: " $name
		perl $api status jobid=$id getResults=no checkEvery=once silent > $lineid.log_status
		set mystatus2 = `grep "ERROR" $lineid.log_status | wc -m`
		if ($mystatus2 > $zero) then
		  echo "Jod $id error-ed."
      rm -f $dir/$name.fasta $lineid*
		  # mv $dir/$name.fasta ${dir}_error/.
		  # mv -f $lineid* ${dir}_error/.
		else
		  set mystatus = `grep "finished. Results available at the following" $lineid.log_status | wc -m`
		  if ($mystatus > $zero) then
  			echo "  Jod $name finished. Getting results, creating job directory and moving there all the relevant files."
        set myurl = `grep "www.compbio.dundee.ac.uk" $lineid.log_status`
        set myhtml = `echo $myurl | sed 's%results.html%concise%'`
        curl -s $myhtml > ${out}/jpred_tmp/$name
        rm -f $dir/$name.fasta $lineid*
        # or mv $lineid* ${dir}_output/.
		  endif
		endif
	end
	#---------------------------------------------------------
  set moreLogs = `ls $dir/*.job_id | wc -l`
end
