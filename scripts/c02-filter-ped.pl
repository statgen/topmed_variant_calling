#!/usr/bin/perl -w

use strict;

## Remove pedigree with
## 0.05 het discordance with family
## 0.01 dup discoraance between dups or
## <95% DP10 covered

#my $ped = "hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.related1.inferred.ped";
my $ped = $ARGV[0];
#my $prefix = "hgdp/merged.autosomes.gtonly.minDP0.hgdp.mendel";
my $prefix = $ARGV[1];

### Read QC file
#open(IN,"manifest/geno/freeze9.manifest.vb_xy.2019_10_25.ALL.dp15_freemix10.tsv") || die "Cannot open file\n";
open(IN, $ARGV[2]) || die "Cannot open manifest file\n";
my $dummy = <IN>;
my %hqc = ();
while(<IN>) {
    #my ($sampleID,$loc,$bucket,$region,$cram,$crai,$prj,$center,$study,$vbdepth,$freemix,$pc1,$pc2,$pc3,$pc4,$cpc1,$cpc2,$cpc3,$cpc4,$fdp1,$fdp5,$fdp10) = split;
    my ($sampleID,$cram,$bcf,$vbdepth,$freemix,$pc1,$pc2,$pc3,$pc4,$cpc1,$cpc2,$cpc3,$cpc4,$fdp1,$fdp5,$fdp10) = split;
    #$hqc{$sampleID} = [$prj,$center,$study,sprintf("%.2lf",$vbdepth),sprintf("%.2g",$freemix),$pc1,$pc2,$pc3,$pc4,sprintf("%.5lf",$fdp10)];
    $hqc{$sampleID} = ["prj","center","study",sprintf("%.2lf",$vbdepth),sprintf("%.2g",$freemix),$pc1,$pc2,$pc3,$pc4,sprintf("%.5lf",$fdp10)];
    #print STDERR "$sampleID\t$prj\n" if ( $sampleID =~ /_CCDG/ );
    #hqc{$sampleID} = sprintf("%.5lf",$fdp10);
}
close IN;

my %hdel = ();

open(IN,"$prefix.ind.dup.conc") || die "Cannot open file\n";
$dummy = <IN>;
while(<IN>) {
    my ($id1,$id2,$total,@ns) = split;
    die $id1 unless (defined($hqc{$id1}));
    die $id2 unless (defined($hqc{$id2}));
    my $tot = 0;
    for(my $i=0; $i < 16; ++$i) { $tot += $ns[$i]; }
    my @a1 = @{$hqc{$id1}};
    my @a2 = @{$hqc{$id2}};
    my $prj = &catcon($a1[0],$a2[0]);
    my $ctr = &catcon($a1[1],$a2[1]);
    my $std = &catcon($a1[2],$a2[2]);
    my $minDP =  &min($a1[3],$a2[3]);
    my $maxMIX = &max($a1[4],$a2[4]);
    my $minFDP10 = &min($a1[9],$a2[9]);
    my $disc = $ns[6]+$ns[7]+$ns[9]+$ns[11]+$ns[13]+$ns[14];
    my $miss = $ns[0]+$ns[1]+$ns[2]+$ns[3]+$ns[4]+$ns[8]+$ns[12];
    my $het1 = $ns[8]+$ns[9]+$ns[10]+$ns[11];
    my $het2 = $ns[2]+$ns[6]+$ns[10]+$ns[14];
    my $hetdisc = ($disc+5e-11)/($tot-$miss-$ns[5]-$ns[15]+1e-10);
    if ( $hetdisc > 0.01 ) {
	$hdel{$id1} = 1;
	$hdel{$id2} = 1;
    }
}
close IN;

### Read Trio concordance. Identify IDs with >0.05 het discordance
open(IN,"$prefix.ind.fam.conc") || die "Cannot open file\n";
$dummy = <IN>;
while(<IN>) {
    my ($dad,$mom,$kid,$total,@ns) = split;
    next if ( ( $dad eq "." ) || ( $mom eq "." ) || ( $kid eq "." ) );
    die $dad unless (defined($hqc{$dad}));
    die $mom unless (defined($hqc{$mom}));
    die $kid unless (defined($hqc{$kid}));        
    my @a1 = @{$hqc{$dad}};
    my @a2 = @{$hqc{$mom}};
    my @a3 = @{$hqc{$kid}};
    my $prj = &catcon($a1[0],$a2[0],$a3[0]);
    my $ctr = &catcon($a1[1],$a2[1],$a3[1]);
    my $std = &catcon($a1[2],$a2[2],$a3[2]);
    my $minDP =  &min($a1[3],$a2[3],$a3[3]);
    my $maxMIX = &max($a1[4],$a2[4],$a3[4]);
    my $minFDP10 = &min($a1[9],$a2[9],$a3[9]);
    my $nmiss = 0;
    my $nconc = 0;
    my $ndisc = 0;
    my @nREFs = (0,0);
    my @nALTs = (0,0);
    my @nHETs = (0,0,0);
    my ($allRef,$allAlt) = 0;
    my $tot = 0;
    for(my $i=0; $i < 64; ++$i) {
	$tot += $ns[$i];
	my $d = int($i / 16);
	my $m = int(($i % 16)/4);
	my $k = $i % 4;
	if ( $d * $m * $k == 0 ) { ++$nmiss; }
	elsif ( $d == 1 ) {
	    if ( $m == 1 ) {
		if ( $k == 1 ) { $nconc += $ns[$i]; $allRef = $ns[$i]; }
		else { $ndisc += $ns[$i]; }
	    }
	    elsif ( $m == 2 ) {
		if ( $k == 3 ) { $ndisc += $ns[$i]; }
		else { $nREFs[$k-1] += $ns[$i]; }
	    }
	    else { ## m == 3
		if ( $k == 2 ) { $nconc += $ns[$i]; }
		else { $ndisc += $ns[$i]; }
	    }
	}
	elsif ( $d == 2 ) {
	    if ( $m == 1 ) {
		if ( $k == 1 ) { $nREFs[0] += $ns[$i]; }
		elsif ( $k == 2 ) { $nREFs[1] += $ns[$i]; }
		else { $ndisc += $ns[$i]; }		
	    }
	    elsif ( $m == 2 ) { ## HETHET
		$nHETs[$k-1] += $ns[$i];
	    }
	    else { ## m == 3
		if ( $k == 1 ) { $ndisc += $ns[$i]; }
		else { $nALTs[$k-2] += $ns[$i]; }
	    }
	}
	else { # d ==3
	    if ( $m == 1 ) {
		if ( $k == 2 ) { $nconc += $ns[$i]; }
		else { $ndisc += $ns[$i]; }
	    }
	    elsif ( $m == 2 ) {
		if ( $k == 1 ) { $ndisc += $ns[$i]; }
		else { $nALTs[$k-2] += $ns[$i]; }		
	    }
	    else {
		if ( $k == 3 ) { $nconc += $ns[$i]; $allAlt = $ns[$i]; }
		else { $ndisc += $ns[$i]; }		
	    }
	}
    }
    my $fmiss = sprintf("%.5lf",$nmiss/$tot);
    my $totdisc = sprintf("%.5lf",($ndisc+5e-11)/($nconc+$ndisc+1e-10));
    my $nrfdisc = sprintf("%.5lf",($ndisc+5e-11)/($nconc+$ndisc-$allRef+1e-10));
    my $hetdisc = sprintf("%.5lf",($ndisc+5e-11)/($nconc+$ndisc-$allRef-$allAlt+1e-10));
    my $refra = sprintf("%.5lf",($nREFs[0]+1e-10)/($nREFs[1]+1e-10));
    my $altra = sprintf("%.5lf",($nALTs[0]+1e-10)/($nALTs[1]+1e-10));
    my $hethet = sprintf("%.5lf",($nHETs[1]+1e-10)/($nHETs[0]+$nHETs[2]+1e-10));
    my $hetra = sprintf("%.5lf",($nHETs[0]+1e-10)/($nHETs[2]+1e-10));
    if ( $hetdisc > 0.05 ) {
	$hdel{$dad} = 1;
	$hdel{$mom} = 1;
	$hdel{$kid} = 1;
    }
}
close IN;

### This is a routine to read a pedigree -- need to make a submodule eventually
my %hfam = ();
my %hdup = ();
my %hsex = ();
open(IN,$ped) || die "Cannot open file\n";
while(<IN>) {
    my ($fam,$ind,$dad,$mom,$sex,$pheno) = split;
    my @inds = split(/,/,$ind);
    $hdup{$inds[0]} = \@inds if ( $#inds > 0 );
    $hsex{$inds[0]} = $sex;
    if ( !defined($hfam{$fam}) ) {
	$hfam{$fam} = [ "", "" ];
    }
    if ( ( $mom eq "0" ) && ( $dad eq "0" ) ) {
	if ( $hfam{$fam}->[$sex-1] ne "" ) {
	    die "Conflicing parents: ".join("\t",$fam,$sex,$hfam{$fam}->[$sex-1],@inds)."\n";
	}
	else {
	    $hfam{$fam}->[$sex-1] = $inds[0];	    
	}
    }
    elsif ( ( $mom eq "0" ) || ( $dad eq "0" ) ) {
	die $_;
    }
    else {
	push(@{$hfam{$fam}},$inds[0]);
    }
}
close IN;


## Scan pedigree and see if anyone disqualfies
foreach my $fam (sort keys %hfam) {
    ## scan the family
    my @ids = @{$hfam{$fam}};
    my $pass = 1;

    ### ADDED: ANNOTATE BY PROJECT
    my %hprjs = ();
    foreach my $id (@ids) {
	next if ( $id eq "" );
	$hprjs{$hqc{$id}->[0]} = 1;
	if ( defined($hdup{$id}) ) {
	    foreach my $id2 (@{$hdup{$id}}) {
		$hprjs{$hqc{$id2}->[0]} = 1;		
	    }
	}
    }
    $fam = "$fam\@".join("_",sort keys %hprjs);
    ### END ADDITION:
    
    foreach my $id (@ids) {
	next if ( $id eq "" );
	die "$id\n" unless ( defined($hqc{$id}) );
	$pass = 0 if ( $hqc{$id}->[9] < 0.95 );
	$pass = 0 if ( defined($hdel{$id}) );
    }
    if ( $pass == 1 ) {
	for(my $i=0; $i < @ids; ++$i) {
	    next if ( $ids[$i] eq "" );
	    my $sex = $hsex{$ids[$i]};
	    my $id = defined($hdup{$ids[$i]}) ? join(",",@{$hdup{$ids[$i]}}) : $ids[$i];
	    if ( $i < 2 ) {
		print join("\t",$fam,$id,0,0,$sex,-9)."\n";
	    }
	    else {
		print join("\t",$fam,$id,$ids[0] ? $ids[0] : 0,$ids[1] ? $ids[1] : 0,$sex,-9)."\n";		
	    }
	}
    }
    else {
	print STDERR "Skipping family $fam..\n";
    }
}

sub min {
    my $ret = 1e100;
    foreach my $f (@_) {
	$ret = $f if ( $ret > $f );
    }
    return $ret;
}


sub max {
    my $ret = -1e100;
    foreach my $f (@_) {
	$ret = $f if ( $ret < $f );
    }
    return $ret;
}

sub catcon {
    my %h = ();
    foreach my $f (@_) {
	$h{$f} = 1;
    }
    return join(",",sort keys %h);
}
