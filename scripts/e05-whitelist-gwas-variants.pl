#!/usr/bin/perl -w

use strict;

my $chr = $ARGV[0];

my %hw = ();
my $nw = 0;
open(IN,"resources/gwascatalog.20200204.uniq.rsid.entries.tsv") || die "Cannot open file\n";
while(<IN>) {
    my ($build,$rsno,@F) = split;
    next unless ( $F[0] eq "chr$chr" );
    my @alts = split(/,/,$F[4]);
    next unless ( $F[7] =~ /;CAF=/ );
    my @cafs = split(/,/,$1) if ( $F[7] =~ /;CAF=([^;]+);/ );
    my ($maxaf,$imax) = (0,0);
    for(my $i=1; $i < @cafs; ++$i) {
	if ( ( $cafs[$i] ne "." ) && ( $cafs[$i] > $maxaf ) ) {
	    $imax = $i;
	    $maxaf = $cafs[$i];
	}
    }
    next if ( $imax == 0 );
    $hw{"$F[1]:$F[3]:$alts[$imax-1]"} = 1;
    ++$nw;
}
close IN;

print STDERR "Finished loading $nw variants to be whitelisted\n";

my $vcf = $ARGV[1]; #"release/sites/nowhite/freeze9.merged.chr$chr.filtered.anno.sites.vcf.gz";
my $outprefix  = $ARGV[2]; #"release/sites/freeze9.merged.chr$chr.filtered.anno.gwas.sites";
my $vcfsummary2 = "$ENV{'EXE_PREFIX'}/apigenome/bin/vcf-summary-v2";
my $ref = "resources/ref/hs38DH.fa";
my $dbsnp = "resources/ref/dbsnp_142.b38.vcf.gz";
my @posVcfs = qw(resources/ref/hapmap_3.3.b38.sites.vcf.gz resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz);
    
open(VCF,"zcat $vcf |") || die "Cannot open file\n";
open(OUT1," | $ENV{'EXE_PREFIX'}/htslib/bgzip -c > $outprefix.vcf.gz") || die "Cannot open file\n";
open(OUT2, "| $vcfsummary2 --ref $ref --db $dbsnp --FNRvcf $posVcfs[0] --chr chr$chr --tabix $ENV{'EXE_PREFIX'}/htslib/tabix --bgzip $ENV{'EXE_PREFIX'}/htslib/bgzip > $outprefix.summary_v2") || die "Cannot open file\n";

my ($ngwasPass,$ngwasKeep,$ngwasSwitch) = (0,0,0);
while(<VCF>) {
    if ( /^#/ ) {
	next if ( /^INFO=<ID=GC/ );
	next if ( /^INFO=<ID=GN/ );	
	next if ( /^INFO=<ID=AF/ );	
	next if ( /^INFO=<ID=HWEAF_P/ );
	next if ( /^INFO=<ID=AVG_IF/ );
	print OUT1 $_;
	if ( /^##INFO=<ID=NM1,/ ) {
	    print OUT1 "##INFO=<ID=GWAS,Number=.,Type=String,Description=\"Whitelisted GWAS catalog variant. Value is the original FILTER column, separated by comma instead of semicolon\">\n";
	}
    }
    else {
	my @F = split(/[\t\r\n]/);
	$F[7] =~ s/;AF=.*;HWEAF_P=[^;]+;/;/;
	$F[7] =~ s/;AVG_IF=[^;]+;/;/;
	my $key = "$F[1]:$F[3]:$F[4]";

	if ( defined($hw{$key}) ) {
	    my $oldFilt = $F[6];
	    $oldFilt =~ s/;/,/g;  ## oldFilt contains the old filters
	    my $newFilt = "PASS"; ## newFilt should be PASS in most cases
	    if ( $oldFilt eq "PASS" ) { ++$ngwasPass; } ## if already PASS, that is fine still pass
	    elsif ( ( $oldFilt =~ /SVM/ ) || ( $oldFilt =~ /CEN/ ) || ( $oldFilt =~ /EXHET/ ) || ( $oldFilt =~ /DISC/ ) || ( $oldFilt =~ /CHRXHET/ ) ) {  ## if failed by existing filters
		++$ngwasKeep; ## keep the current filter
		$newFilt = $F[6];
	    }
	    else { ## must be only MIS2,DUP2,TRI2
		++$ngwasSwitch;
	    }
	    $F[7] .= ";GWAS=$oldFilt";
	    $F[6] = $newFilt; #"PASS";
	    print STDERR "$F[0]:$F[1]:$F[3]:$F[4] $ngwasPass $ngwasKeep $ngwasSwitch\n" if ( rand() < 0.01 );
	}
	print OUT1 join("\t",@F)."\n";
	print OUT2 join("\t",@F)."\n";
    }
}
close OUT1;
close OUT2;
close VCF;

print STDERR "Finished $ngwasPass $ngwasKeep $ngwasSwitch\n";

print `$ENV{'EXE_PREFIX'}/htslib/tabix -f -pvcf $outprefix.vcf.gz`;
