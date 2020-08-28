#!/usr/bin/perl -w
  
use strict;

#my $chr = $ARGV[0];

my $milkvcf = $ARGV[0]; #"fixed0/milk.filt/milk.chr$chr.merged.sites.vcf.gz";
my $fmisbcf = $ARGV[1]; #"fixed0/sites.update_info/merged.chr$chr.gtonly.minDP10.update_info.sites.bcf";
#my $out = "/dev/stdout" #"analysis/filt/ld/frz9/frz9.milk_nold.fmis10.chr$chr.vcf.gz";

open(BCF,"$ENV{'EXE_PREFIX'}/bcftools/bcftools view -H $fmisbcf |") || die "Cannot open file\n";

open(IN,"zcat $milkvcf |") || die "Cannot open file\n";
open(OUT,"| $ENV{'EXE_PREFIX'}/htslib/bgzip -c") || die "Cannot open file\n";
while(<IN>) {
    print STDERR "Processing $. lines..\n" if ( $. % 1000000 == 0 );
    if ( /^#/ ) {
        print OUT $_;
        if ( /ID=TRIO_CONC_THRES/ ) {
            print OUT "##INFO=<ID=FMIS10,Number=1,Type=Float,Description=\"Fraction of missing genotype at depth 10\">\n";
        }
    }
    else {
        my @F = split;
        my @B = split(/[\t\r\n ]+/,<BCF>);
        next unless ( ( $F[1] eq $B[1] ) || ( $F[3] eq $B[3] ) || ( $F[4] eq $B[4] ) );
        my $fmis = $1 if ( $B[7] =~ /;FMIS=(\S+)/ );
        $F[7] =~ s/;MILK_LRE=/;FMIS10=$fmis;MILK_LRE=/;
        print OUT join("\t",@F)."\n";
    }
}
close OUT;
close IN;
close BCF;

#print `tabix -f -pvcf $out`;
