#!/usr/bin/perl -w

## Cloudify is a script that helps users to run jobs on the cloud

## this time, put all job status into a GCS bucket
## the object name should be the same for the same task
## (i.e. only dependent on parameters)
## the object should be constantly updated to record
## [jobname, machine name, last timestamp checked, last timestamp the current job started, startup script, the sequence of current job]
## Before creating the machine, the object must be created

use warnings;
use POSIX qw(strftime);
use strict;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(forkExecWait);
use wGetOptions qw(wpod2usage wGetOptions);

my $outDir         = ""; ## local directory to store output files 
my $jobName        = ""; ## name of the job name. $jobname.sh will be the script name, machine name will be $jobname
my $imageName      = ""; ## "ubuntu-20170426-vt";
my @zones          = (); #"us-central1-".$zones[int(rand($#zones+1))];
my $preemptible    = "";
my $opts           = ""; ## "--custom-cpu 1 --custom-memory 6656MiB";
my $disk           = 200;
my @cmds           = ();
my @mounts         = ();
my $gsdir          = "";  ## GCS bucket directory to store metadata and scripts
my $keepInstance   = "";
my $cwd            = "";
my $checkpointSecs = 0;
my $maxHours       = 0;
my $sleepSecs      = 0;
my @targetFiles    = ();
my %lists          = ();

open(CMD,$ARGV[0]) || die "Cannot open $ARGV[0]\n";
while(<CMD>) {
    if ( /^#/ ) {
	## this part would 'OVERRIDE' existing default value
	if ( /^#+\s*(\S+)\s*[:=]\s*(\S.*\S|\S)\s*$/ ) {
	    my ($key,$value) = ($1,$2);
	    chomp $value;
	    if ( $key eq "list" ) {
		if ( $value =~ /^(\S+)\s*[:=]\s*(\S.*\S|\S)$/ ) {
		    my ($alias,$file) = ($1,$2);
		    unless ( defined($file) && $file ) {
			die "Cannot recognize list value $value\n";
		    }
		    my @alist = ();
		    my $ncol = 0;
		    open(IN,$file) || die "Cannot open file $file\n";
		    while(<IN>) {
			next if ( /^#/ );
			my @F = split;
			push(@alist,\@F);
			if ( $ncol == 0 ) { $ncol = $#F+1; }
			elsif ( $ncol != $#F+1 ) {
			    die "The expected number of columns is $ncol, but ".($#F+1)." is observed\n";
			}
		    }
		    close IN;
		    $lists{$alias} = \@alist;
		}
	    }
	    elsif ( $key eq "target" ) {
		if ( $#targetFiles >=0 ) { print STDERR "Overriding --target argument to $value...\n"; }
		unless ( $value =~ /^gs:\/\// ) {
		    die "Cannot recognize target value $value. Expected to start with gs://\n";
		}
		@targetFiles = split(/\s+/,$value);
		#push(@targetFiles, $value);
		#$targetFile = $value;
	    }	    
	    elsif ( $key eq "out" ) {
		if ( $outDir ) { print STDERR "Overriding --out argument to $value...\n"; }
		$outDir = $value;		    
	    }
	    elsif ( $key eq "name" ) {
		if ( $jobName ) { print STDERR "Overriding --name argument to $value...\n"; }
		$jobName = $value;
	    }
	    elsif ( $key eq "image" ) {
		if ( $imageName ) { print STDERR "Overriding --image argument to $value...\n"; }
		$imageName = $value;		    
	    }
	    elsif ( $key eq "gs" ) {
		if ( $gsdir ) { print STDERR "Overriding --gs argument to $value...\n"; }
		$gsdir = $value;		    		    
	    }
	    elsif ( $key eq "disk" ) {
		if ( $disk ) { print STDERR "Overriding --disk argument to $value...\n"; }
		unless ( $disk =~ /^\d+$/ ) { die "disk parameter must be integer (in GB) in $cmds[0]\n"; }
		$disk = $value;		    		    
	    }
	    elsif ( $key eq "zone" ) {
		if ( $#zones >= 0 ) { print STDERR "Overriding --zone argument to $value...\n"; }
		@zones = split(/\s+/,$value);
		#$zone = $value;		    		    		    
	    }
	    elsif ( $key eq "cwd" ) {
		if ( $cwd ) { print STDERR "Overriding --name argument to $value...\n"; }
		$cwd = $value;		    		    		    		    
	    }
	    elsif ( $key eq "mount" ) {
		my @values = split(/\s+/,$value);
		if ( $#mounts >= 0 ) { print STDERR "Overriding --mount argument to @values (array) ...\n"; }
		@mounts = @values;
	    }
	    elsif ( $key eq "opts" ) {
		if ( $opts ) { print STDERR "Overriding --name argument to $value...\n"; }
		$opts = $value;		    		    		    		    
	    }
	    elsif ( $key eq "preemptible" ) {
		if ( ( $value eq "true" ) || ( $value eq "TRUE" ) ) { $value = 1; }
		elsif ( ( $value eq "false" ) || ( $value eq "FALSE" ) ) { $value = 0; }
		elsif ( ( $value ne 0 ) && ( $value ne 1 ) ) {
		    die "--preemptible argument must be either true or false, but the current value is '$value'\n";
		}
		if ( $preemptible ) { print STDERR "Overriding --preemtible argument to $value...\n"; }
		$preemptible = $value;		    		    		    		    
	    }
	    elsif ( $key eq "checkpoint-secs" ) {
		if ( $checkpointSecs ) { print STDERR "Overriding --checkpoint-secs argument to $value...\n"; }
		unless ( $value =~ /^\d+$/ ) { die "--checkpoint-secs parameter must be integer (in GB) in $cmds[0]\n"; }
		$checkpointSecs = $value;		    		    		    
	    }
	    elsif ( $key eq "max-hours" ) {
		if ( $maxHours ) { print STDERR "Overriding --max-hours argument to $value...\n"; }
		unless ( $value =~ /^\d+$/ ) { die "--max-hours parameter must be integer in $cmds[0]\n"; }
		$maxHours = $value;		    		    		    
	    }
	    elsif ( $key eq "sleep-secs" ) {
		if ( $sleepSecs ) { print STDERR "Overriding --sleep-secs argument to $value...\n"; }
		unless ( $value =~ /^\d+$/ ) { die "--sleep-secs parameter must be integer in $cmds[0]\n"; }
		$sleepSecs = $value;		    		    		    
	    }		    
	    elsif ( $key eq "keep" ) {
		if ( ( $value eq "true" ) || ( $value eq "TRUE" ) ) { $value = 1; }
		elsif ( ( $value eq "false" ) || ( $value eq "FALSE" ) ) { $value = 0; }
		elsif ( ( $value != 0 ) && ( $value != 1 ) ) {
		    die "--keep argument must be either true or false in $cmds[0]\n";
		}
		if ( $keepInstance ) { print STDERR "Overriding --keep argument to $value...\n"; }
		$keepInstance = $value;		    		    		    		    		    
	    }
	    else {	
		print STDERR "Cannot recognize --$key argument to replace into the value $value. Skipping...\n";
	    }
	}
	else {
	    print STDERR "Cannot parse the header line $_, skipping..\n";
	}
    }
    elsif ( /^\s*$/ ) { next; }
    else {
	chomp;
	s/^\s+//g;
	s/\s+$//g;
	push(@cmds,$_);
    }
}
close CMD;

print STDERR "Successfully read ".($#cmds+1)." commands from $ARGV[0]\n";

unless ( ( $jobName ) && ( $imageName ) && ( $#cmds >= 0 ) && ( $disk > 0 ) && ( $gsdir =~ /^gs:\/\// ) && ( $#zones >= 0 ) && ( $imageName ) ) {
    print STDERR "Missing required options (and make sure that $gsdir is a valid GCS URL\n";
    wpod2usage(2);
}

print STDERR "------------------------------------------------\n";
print STDERR "Parameters in effect:\n";
print STDERR "------------------------------------------------\n";
print STDERR "name:\t$jobName\n";
print STDERR "image:\t$imageName\n";
print STDERR "cmd:\t(Total ".($#cmds+1).")\n";
for(my $i=0; $i < @cmds; ++$i) {
    print STDERR "  ".($i+1).": $cmds[$i]\n";
}
print STDERR "gs:\t$gsdir\n";
print STDERR "disk:\t$disk GB\n";
print STDERR "zones:\t@zones\n";
print STDERR "cwd:\t$cwd\n";
print STDERR "mount:\t(Total ".($#mounts+1).")\n";
for(my $i=0; $i < @mounts; ++$i) {
    print STDERR "  ".($i+1).": $mounts[$i]\n";
}
print STDERR "opts:\t$opts\n";
print STDERR "checkpoint-secs:\t$checkpointSecs seconds\n";
print STDERR "max-hours:\t$maxHours hours\n";
print STDERR "preemptible:\t".($preemptible ? "true" : "false")."\n";
print STDERR "keep:\t".($keepInstance ? "true" : "false")."\n";
print STDERR "------------------------------------------------\n";

my @listkeys = sort keys %lists;
my @listvals = ();
foreach my $key (@listkeys) {
    push(@listvals,$lists{$key});
}
my @listprods = &array_prod(@listvals);

print STDERR "-------------------------------------------------\n";
print STDERR "List items to iterate:\n";
print STDERR "-------------------------------------------------\n";
foreach my $listkey (@listkeys) {
    my $r = $lists{$listkey};
    my $nitems = $#{$r}+1;
    if ( $nitems == 0 ) {
	die "Error: No item observed for list key $listkey\n";
    }
    my $ncol = $#{$r->[0]}+1;
    print STDERR "$listkey:\t$nitems items, $ncol columns\n";
}
print STDERR "-------------------------------------------------\n";
print STDERR "Total: :".($#listprods+1)." items\n";
print STDERR "-------------------------------------------------\n";

my $cmdstr = "(\"".join("\",\"",@cmds)."\")";
my $mountstr = ($#mounts < 0) ? "()" : "(\"".join("\",\"",@mounts)."\")";

&forkExecWait("mkdir --p $outDir") unless ( -e $outDir );

my $msg =<<"END_MESSAGE1";
#!/usr/bin/perl -w

use warnings;
use POSIX qw(strftime :sys_wait_h);
use strict;

my \$cwd = "$cwd";          # working directory
my \@mounts = $mountstr;
my \@cmds = $cmdstr;
my \$gsdir = "$gsdir";
my \$keep = "$keepInstance";        # delete the GCE VM instance upon successful execution
my \$checkpointSecs = $checkpointSecs;
my \$maxHours = $maxHours;
END_MESSAGE1

$msg .=<<'END_MESSAGE2';
## We assume that the machine was already created.
my $vmName = `curl "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google"`;
unless ( $vmName ) {
    print STDERR "Cannot determine machine name \n";
    unless ( $keep ) { &shutdownVM(); }
}

$gsdir =~ s/\/$//g;
my $gsLogFile  = "$gsdir/$vmName.log";
my $gsErrFile  = "$gsdir/$vmName.err";

my $vmZone = $1 if (`curl "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google"` =~ /zones\/(\S+)$/ );

my $pid = $$;
my $ts = time();

unless ( $vmZone ) {
    print STDERR "ERROR: Unknown Zone from $vmName\n";
    &transferLog($gsErrFile);
    unless ( $keep ) { &shutdownVM(); }
}

## determine the disk names to be mounted
my @diskdevs = qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);
my $ilast = -1;
for(my $i=0; $i < @diskdevs; ++$i) {
    my $l = $diskdevs[$i];
    if ( -e "/dev/sd$l" ) { $ilast = $i; }
    else { last; }
}

### attach and mout each disks
for(my $i=0; $i < @mounts; ++$i) {
    my ($diskname,$path,$mode) = split(/,/,$mounts[$i]);
    $mode = "rw" unless ( $mode );
    
    my $cmd = "gcloud compute instances attach-disk $vmName --disk $diskname --zone $vmZone --mode $mode";
    my $ret = &forkExecWait($cmd);
    if ( $ret ) {
      print STDERR "ERROR: Could not attach disk $diskname\n";
      &transferLog($gsErrFile);
      unless ( $keep ) { &shutdownVM(); }

    }
    print STDERR "Successfully attached disk $diskname\n";
    
    my $l = $diskdevs[$ilast+$i+1];
    $cmd = "sudo mkdir --p $path && sudo mount -o discard,defaults /dev/sd$l $path";
    $ret = &forkExecWait($cmd);
    if ( $ret ) { 
      print STDERR "ERROR: Could not mount disk $diskname at /dev/sd$l to $path\n";
      &transferLog($gsErrFile);
      unless ( $keep ) { &shutdownVM(); }
    }

    print STDERR "Successfully mounted disk $diskname at /dev/sd$l to $path\n";
}

if ( $cwd ) { 
    if ( chdir($cwd) ) {
        print STDERR "Successfully changed the working directory to $cwd\n";
    }
    else {
        die "ERROR: Could not change the working directory to $cwd\n";
    }
}


if ( $#mounts >= 0 ) {
  print STDERR "Successfully mounted the disks\n";
}

for(my $i=0; $i < @cmds; ++$i) {
    print STDERR "Running command '$cmds[$i]' on target...\n";
    my $cmd = $cmds[$i];
    if ( $checkpointSecs > 0 ) {
      my $ret = &forkExecCheckpoint($cmd,$checkpointSecs,$maxHours*3600);
      if ( $ret ) { 
         print STDERR "ERROR: Could not successfully finish the command '$cmd', returning exit code $ret\n";
         &transferLog($gsErrFile);
         unless ( $keep ) { &shutdownVM(); }
      }
    }
    else {
      my $ret = &forkExecWait($cmd);
      if ( $ret ) { 
         print STDERR "ERROR: Could not successfully finish the command '$cmd', returning exit code $ret\n";
         &transferLog($gsErrFile);
         unless ( $keep ) { &shutdownVM(); }
      }
    }
    print STDERR "Successfully finished the step$i command '$cmd'\n";
}

$ts = time(); 

&transferLog($gsLogFile);
unless ( $keep ) { &shutdownVM(); }

sub transferLog {
    my $logfile = shift;
    my $cmd = "sudo gsutil cp /var/log/startupscript.log $logfile";
    my $ret = &forkExecWait($cmd);
    if ( $ret ) { die "ERROR: Could not copy the log file of compute instance $vmName in $vmZone\n"; }
    else { print STDERR "Successfully copied the syslog before deleting the machine\n"; }
}

sub shutdownVM {
    my $cmd = "gcloud compute instances delete $vmName --zone $vmZone --quiet";
    my $ret = &forkExecWait($cmd);
    if ( $ret ) { die "ERROR: Could not delete compute instance $vmZone\n"; }
    else { print STDERR "Successfully shutting down the machine $vmName in $vmZone\n"; }
}

sub forkExecWait {
    my $cmd = shift;
    my $msg = shift;
    if ( defined($msg) ) {
	print "$msg...";
    }
    else {
	print "forkExecWait(): $cmd ...";
    }
    my $kidpid;
    if ( !defined($kidpid = fork()) ) {
	die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
	exec($cmd);
	die "Cannot exec $cmd: $!";
    }
    else {
	waitpid($kidpid,0);
    }

    print ".. finished!\n";
    return ($?>>8);
}

sub forkExecCheckpoint {
    my ($cmd,$period,$timeout) = @_;

    print "forkExecCheckpoint():\ncmd: $cmd\nperiod: $period\ntimeout: $timeout\n";
    my $kidpid;
    my $startTime = time();
    if ( !defined($kidpid = fork()) ) {
        die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
        exec($cmd);
        die "Cannot exec $cmd: $!";
    }
    else {
	for(my $t=0;$t < $timeout;) {
	    my $ret = waitpid($kidpid,WNOHANG);
	    if ( $ret == -1 ) {
		my $errcode = ($? >> 8);
		print STDERR "... error occurred with code $errcode\n";
		return $errcode;
	    }
	    elsif ( $ret ) {
		my $errcode = ($? >> 8);
		print STDERR "... finished with exit code $errcode\n";
		return $errcode;		
	    }
  
            print "ts=$t, kidpid=$kidpid, $timeout = $timeout, ret=$ret\n";

            if ( $t == 0 ) { sleep 5; $t += 5; }
            elsif ( $t == 5 ) { sleep 60; $t += 60; }
	    else { sleep $period; $t += $period; }
	}
    }
    print ".. finished!\n";
    return ($?>>8);
}
END_MESSAGE2

for(my $i=0; $i < @listprods; ++$i) {
    ## substitute variables
    my $myMsg = subst_keywords(\@listkeys,$listprods[$i],$msg);
    open(OUT,">$outDir/$jobName.$i.pl") || die "Cannot open $outDir/$jobName.$i.pl for writing\n";
    print OUT $myMsg;
    close OUT;
}

$gsdir =~ s/\/$//g;

open(MAK,">$outDir/$jobName.mk") || die "Cannot open file\n";
print MAK ".DELTE_ON_ERROR:\n\n";
print MAK "all:";
for(my $i=0; $i < @listprods; ++$i) {
    print MAK " $outDir/$jobName.$i.OK";
}
print MAK "\n\n";
for(my $i=0; $i < @listprods; ++$i) {
    print MAK "$outDir/$jobName.$i.OK:\n";
    ## check if targetFile exists
    my @myTargets = ();
    for(my $j=0; $j < @targetFiles; ++$j) {
	push(@myTargets, subst_keywords(\@listkeys,$listprods[$i],$targetFiles[$j]));
    }
    my $machineName = lc("$jobName-$i");
    my $ts = time();
    my $zone = $zones[$i % ($#zones+1)];
    print MAK "\t(gsutil ls @myTargets > /dev/null 2> /dev/null && touch $outDir/$jobName.$i.OK && echo 'target @myTargets already exists') || (gcloud compute instances describe --zone $zone $machineName > /dev/null 2> /dev/null && echo '$machineName is already running') || (gsutil cp $outDir/$jobName.$i.pl $gsdir/$jobName.$i.pl && gcloud compute instances create --zone $zone --scopes cloud-platform --image $imageName $opts --boot-disk-size $disk ".($preemptible ? "--preemptible " : "")." --metadata startup-script-url=$gsdir/$jobName.$i.pl $machineName && sleep $sleepSecs)\n";
    print MAK "\n";
}
close MAK;

print "------------------------------------------------------------------\n";
print "To test, run 'make -f $outDir/$jobName.mk' to test a few cases\n";
print "To run all the jobs, run 'make -f $outDir/$jobName.mk -k -j 10'\n";
print "------------------------------------------------------------------\n";

sub array_prod {
    my ($curlist,@rest) = @_;
    my @curprod = ();
    if ( $#rest >= 0 ) {
	@curprod = array_prod(@rest);  ## curprod = ( [ ["c",3] ], [ ["d",4] ] )
    }
    else {
	@curprod = ( [] );
    }

    my @newprod = ();
    foreach my $item (@{$curlist}) { ## $item = ["a",1]
	foreach my $item2 (@curprod) { ## #item2 = [ ["c", 3] ]
	    push(@newprod,[$item,@{$item2}]);
	}
    }
    return @newprod;
}

## substitute keywords from string
## 
sub subst_keywords {
    my ($rKey,$rProd,$msg) = @_;
    my @keys = @{$rKey};
    my @prod = @{$rProd};
    if ( $#keys != $#prod ) {
	die "Keys and Prod does not have the same length : $#keys +1 != $#prod +1\n";
    }
    for(my $i=0; $i < @keys; ++$i) {
	my @vals = @{$prod[$i]};
	for(my $j=0; $j < @vals; ++$j) {
	    my $k = $j+1;
	    # $KEYWORD$1$
	    #print STDERR "\$$keys[$i]\$$k\$\t$vals[$j]\n";
	    $msg =~ s/\$$keys[$i]\$$k\$/$vals[$j]/g;
	}
    }
    return $msg;
}
