#!/usr/bin/perl -w

use strict;

my $chr = $ARGV[2];

my $milk = $ARGV[1]; #"fixed0/milk.filt/milk.chr$chr.merged.sites.vcf.gz";
my $svm  = $ARGV[0]; #"analysis/filt/ld/frz9/svm.hm3ld.fmis10/frz9.milk_svm.hm3ld_fmis10.chr$chr.sites.vcf.gz";
my $outprefix  = $ARGV[3]; #"fixed0/svm.ld.fmis/frz9.milk_svm.release.chr$chr.sites";
my $vcfsummary2 = "$ENV{'EXE_PREFIX'}/apigenome/bin/vcf-summary-v2";
my $ref = "resources/ref/hs38DH.fa";
my $dbsnp = "resources/ref/dbsnp_142.b38.vcf.gz";
my @posVcfs = qw(resources/ref/hapmap_3.3.b38.sites.vcf.gz resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz);
my ($xLabel,$xStart,$xStop) = ("chrX", 2781479, 155701383);
    
open(SVM,"zcat $svm|") || die "Cannot open file\n";
open(MILK,"zcat $milk | grep -v ^# |") || die "Cannot open file\n";
open(OUT1," | $ENV{'EXE_PREFIX'}/htslib/bgzip -c > $outprefix.vcf.gz") || die "Cannot open file\n";
open(OUT2, "| $vcfsummary2 --ref $ref --db $dbsnp --FNRvcf $posVcfs[0] --chr $chr --tabix $ENV{'EXE_PREFIX'}/htslib/tabix --bgzip $ENV{'EXE_PREFIX'}/htslib/bgzip > $outprefix.summary_v2") || die "Cannot open file\n";
while(<SVM>) {
    if ( /^#/ ) {
	next if ( /^##FILTER=<ID=DISC,/ );	
	print OUT1 $_;
        #if ( /^##INFO=<ID=LOG_HM3_DIST_08/ ) {
	if ( /^##INFO=<ID=FMIS/ ) {
	    print OUT1 "##INFO=<ID=DUP_NH_ALL,Number=1,Type=Integer,Description=\"Number of duplicate pairs where both genotype depths are 10+ and not identical homozygotes\">\n";
	    print OUT1 "##INFO=<ID=DUP_NH_DIS,Number=1,Type=Integer,Description=\"Number of discordant duplicate genotype pairs where both depths are 10+ and not identical homozygotes\">\n";
	    print OUT1 "##INFO=<ID=TRI_NH_ALL,Number=1,Type=Integer,Description=\"Number of trios where all genotype depths are 10+, mendelian-informative, and not identical homozygotes\">\n";
	    print OUT1 "##INFO=<ID=TRI_NH_DIS,Number=1,Type=Integer,Description=\"Number of mendelian-discordant trios where all genotype depths are 10+, mendelian-informative, and not identical homozygotes\">\n";
	    print OUT1 "##FILTER=<ID=DUP2,Description=\"Duplicate genotype discordance is greater than 2%, with at least two discordances\">\n";
	    print OUT1 "##FILTER=<ID=TRI2,Description=\"Mendelian genotype discordance is greater than 2%, with at least two discordances\">\n";
	    print OUT1 "##FILTER=<ID=MIS2,Description=\"Genotype missing rate at depth 10 is greater than 2%\">\n";	    
	}
    }
    else {
	my @F = split(/[\t\r\n]/);
	my @M = split(/[\t\r\n]/,<MILK>);
	die unless ( $F[1] == $M[1] );

	print STDERR "Processing $F[0]:$F[1]\n" if ( $. % 1000000 == 0 );

	my ($dup,$trio) = ($1,$2) if ( $M[7] =~ /;DUP_CONC_THRES=([^;]+);.*;TRIO_CONC_THRES=([^;]+)/ );
	my $fmis = $1 if ( $F[7] =~ /;FMIS10=([^;]+)/ );
	my @dups = split(/,/,$dup);
	my @trios = split(/,/,$trio);

        ## calculate dup discordance
	my $duphet  = $dups[4];
	my $dupdisc = $dups[1]+$dups[2]+$dups[3]+$dups[5]+$dups[6]+$dups[7];
	my $dupnhom = $duphet + $dupdisc;
	my $dupnref  = $dupnhom + $dups[8];
	my $dupall   = $dupnref + $dups[0];
	my $dupFilt = ( ($dupdisc > 1 ) && ( $dupdisc > 0.02*$dupnhom) ) ? 1 : 0;

	## 0,1,2    : 0,0,0 C  0,0,1 D  0,0,2 D
	## 3,4,5    : 0,1,0 A  0,1,1 A  0,1,2 D
	## 6,7,8    : 0,2,0 D  0,2,1 C  0,2,2 D
	## 9,10,11  : 1,0,0 A  1,0,1 A  1,0,2 D
	## 12,13,14 : 1,1,0 A  1,1,1 A  1,1,2 A
	## 15,16,17 : 1,2,0 D  1,2,1 A  1,2,2 A
	## 18,19,20 : 2,0,0 D  2,0,1 C  2,0,2 D
	## 21,22,23 : 2,1,0 D  2,1,1 A  2,1,2 A
	## 24,25,26 : 2,2,0 D  2,2,1 D  2,2,2 C 
	my @idxD = (1,2,5,6,8,11,15,18,20,21,24,25);
	my @idxC = (0,7,19,26);

        ## for chrX
        ## 0,1,2    : 0,0,0 C  0,0,1 D  0,0,2 D
        ## 3,4,5    : 0,1,0 A  0,1,1 A  0,1,2 A
        ## 6,7,8    : 0,2,0 A  0,2,1 C  0,2,2 C
        ## 9,10,11  : 1,0,0 A  1,0,1 A  1,0,2 A
        ## 12,13,14 : 1,1,0 A  1,1,1 A  1,1,2 A
        ## 15,16,17 : 1,2,0 A  1,2,1 A  1,2,2 A
        ## 18,19,20 : 2,0,0 C  2,0,1 C  2,0,2 A
        ## 21,22,23 : 2,1,0 A  2,1,1 A  2,1,2 A
        ## 24,25,26 : 2,2,0 D  2,2,1 D  2,2,2 C
        my @idxDX = (1,2,24,25);
        my @idxCX = (0,7,8,18,19,26);

        if ( ( $xLabel eq $chr ) && ( $F[1] >= $xStart ) && ( $F[1] <= $xStop ) ) {
            @idxD = @idxDX;
            @idxC = @idxCX;
        }
        print STDERR "$F[1] @idxD $xLabel $chr\n" if ( $F[1] % 1000000 == 0 );

	my ($trioconc,$triodisc) = (0,0);
	foreach my $i (@idxD) { $triodisc += $trios[$i]; }
	foreach my $i (@idxC) { $trioconc += $trios[$i]; }
	
	my $trioall = $trioconc + $triodisc;
	my $trionref = $trioall - $trios[0];
	my $trionhom = $trionref - $trios[26];
	my $triFilt = ( ($triodisc > 1 ) && ( $triodisc > 0.02*$trionhom) ) ? 1 : 0;
	my $misFilt = ( $fmis > 0.02 ) ? 1 : 0;
	
	my @filts = split(/;/,$F[6]);
	my @newfilts = ();
	foreach my $f (@filts) {
	    push(@newfilts,$f) if ( ( $f ne "PASS" ) && ( $f ne "DISC" ) );
	}
	push(@newfilts,"DUP2") if ( $dupFilt == 1 );
	push(@newfilts,"TRI2") if ( $triFilt == 1 );	
	push(@newfilts,"MIS2") if ( $misFilt == 1 );
	push(@newfilts,"PASS") if ( $#newfilts < 0 );

	$F[7] =~ s/;SVM=/;DUP_NH_ALL=$dupnhom;DUP_NH_DIS=$dupdisc;TRI_NH_ALL=$trionhom;TRI_NH_DIS=$triodisc;SVM=/;
	print OUT1 join("\t",@F[0..5],join(";",@newfilts),$F[7])."\n";
	print OUT2 join("\t",@F[0..5],join(";",@newfilts),$F[7])."\n";	
    }
}
close OUT1;
close OUT2;
close MILK;
close SVM;

print `$ENV{'EXE_PREFIX'}/htslib/tabix -f -pvcf $outprefix.vcf.gz`;
