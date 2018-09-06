#!/usr/bin/perl -w                                                                                                                                                          
use warnings;
use POSIX qw(strftime);
use strict;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(forkExecWait makeMake);
use wGetOptions qw(wpod2usage wGetOptions);

my @szchrs = (248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569);
my %hszchrs = ();
my @achrs = (1..22,"X","Y","M");

for(my $i=0; $i < @szchrs; ++$i) {
    $hszchrs{"chr$achrs[$i]"} = $szchrs[$i];
}

my @chrs = @ARGV;
my $prefix = "tmp/vt.site";
my $imageName = "ubuntu-20170610a-vt";
my $zone = "us-central1-b";
my $opts = "--custom-cpu 1 --custom-memory 6656MiB";
my $disk = 200;

my $out = "out.20170521";

foreach my $chr (@chrs) {
    #my $cmd = "gsutil mkdir gs://topmed-vt/paste/$chr";
    #&forkExecWait($cmd);
    my $sz = $hszchrs{$chr};

    for(my $beg=1; $beg < $sz; $beg += 1000000) {
	my $end = $beg + 999999 < $sz ? $beg + 999999 : $sz;
	my $lchr = lc($chr);	
	my $gsout = `gsutil ls -l gs://topmed-vt/paste/$out/$chr/$chr\_$beg\_$end.genotypes.bcf.csi`;
	if ( $? >> 8 == 0 ) {
	    $gsout = `gsutil ls -l gs://topmed-vt/paste/$out/$chr/$chr\_$beg\_$end.unfiltered.sites.bcf`;
	    if ( $? >> 8 != 0 ) { ## the site does not exist
		open(CMD,">$prefix.$chr.$beg.$end.cmd") || die "Cannot open file\n";
		print CMD "# script: $prefix.$chr.$beg.$end.sh\n";
		print CMD "# image: $imageName\n";
		#			print CMD "# cwd: /mnt/indices/crams\n";
		print CMD "# log: gs://topmed-vt/paste/$out/logs\n";
		print CMD "# disk: 100\n";
		print CMD "# zone: us-central1-b\n";
		print CMD "# name: vt-site-c-$lchr-$beg-$end\n";
		#			print CMD "# mount: topmed-cram-indices,/mnt/indices,ro topmed-sites-scripts,/mnt/sites,ro\n";
		print CMD "# opts: --custom-cpu 1 --custom-memory 6656MiB\n";
		print CMD "# prefix: $prefix\n";
		print CMD "# check-every-second: 1200\n";
		print CMD "# check-for-hours: 120\n";
		print CMD "# preemptible: true\n";
		print CMD "# keep: false\n";
		print CMD "chown -R hmkang:hmkang /home/hmkang\n";
		print CMD "sudo -u hmkang gcsfuse --implicit-dirs topmed-vt /mnt/topmed-bcf\n";
		print CMD "sudo -u hmkang bcftools view -G /mnt/topmed-bcf/paste/$out/$chr/$chr\_$beg\_$end.genotypes.bcf -Ob -o /home/hmkang/$chr\_$beg\_$end.unfiltered.sites.bcf\n";
		print CMD "sudo -u hmkang gsutil cp /home/hmkang/$chr\_$beg\_$end.unfiltered.sites.bcf gs://topmed-vt/paste/$out/$chr/\n";		
		close CMD;
			
		my $cmd = "nohup time perl cloudify2.pl --cmd $prefix.$chr.$beg.$end.cmd > $prefix.$chr.$beg.$end.cmd.nohup &";
		print "$cmd\n";
		&forkExecWait($cmd);
	    }
	}
	else {
	    print STDERR "Skipping unfinished $prefix.$chr.$beg.$end.cmf";
	#		$gsout = `gcloud compute instances describe vt-geno-e-$lchr-$beg-$end --zone us-central1-b 2> /dev/null`;
	#		if ( $? >> 8 == 0 ) {			
	#		    print STDERR "Skipping running $prefix.$chr.$beg.$end.cmd..\n";
	#		}
	#		else {
	#		    print STDERR "ERROR: Neither vt-geno-d-$lchr-$beg-$end nor vt-geno-d-$lchr-$beg-$end is running for unfinished chunk\n"; 
	#		}
	#	    }
	#	}
	}
        #}
	#else {
	#    print STDERR "Cannot find BCF sites for $prefix.$chr.$sbeg.$send.cmd. Exit code is $?\n"; #$gsout\n";
	#}
	#sleep 1;
    }
}
#&forkExecWait($cmd);

#}
