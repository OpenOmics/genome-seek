#!/usr/bin/perl

use List::Util qw(min max);
use Getopt::Long;
use strict;

my $svtools = "/opt/conda/envs/python2/bin/svtools";
my $bedtools = "/usr/local/bin/bedtools";
my $minimap = "/usr/local/bin/minimap2";
my $tee = "/usr/bin/tee";
my $refseq = "/gscuser/dspencer/refdata/hg38/all_sequences.fa";
my $tmp = "/tmp";

my $slop = 200;
my $masklevel = 0.8;
my $PR = 2;
my $SR = 2;
my $minSVabund = 5.0;
my $delcovratio = 0.8;
my $dupcovratio = 1.3;
    
my $minMQ = 0;
my $fracIdent = 0.95;

my $knownfile = '';

my $usage = <<END;
$0 -k <knowntrans bedpe> <gzipped vcf> <out vcf name>

  options:
    -r refseq
    
    [ manta hit filtering ] 
    -m masklevel
    -pr PR reads
    -sr SR reads
    -a min VAF/SVabund
    -l slop length for known translocations
    -q min map quality of contig hits
    -i min fraction identity for contig hits

    [ file/tool paths ]
    -b bedtools path
    -s svtools path
    -p minimap path
    -t tmpdir 

END

GetOptions("pr=i" => \$PR,
	   "sr=i" => \$SR,
	   "a=f" => \$minSVabund,
	   "b=s" => \$bedtools,
	   "s=s" => \$svtools,
	   "p=s" => \$minimap,
	   "t=s" => \$tmp,
	   "m=f" => \$masklevel,
           "r=s" => \$refseq,
	   "l=i" => \$slop,
	   "q=i" => \$minMQ,
	   "i=f" => \$fracIdent,
	   "k=s" => \$knownfile);


die "$usage\n" if ! -e $knownfile || ! -e $ARGV[0] || !$ARGV[1] || ! -e $svtools || ! -e $bedtools || ! -e $minimap || ! -e $refseq || $ARGV[0] !~ /vcf.gz$/;

my $out = $ARGV[1];

#open manta VCF
open(VCF,"/bin/zcat $ARGV[0] |") || die "Cant open manta file $ARGV[0]";

# create a fastq file
open(FQ,">$tmp/input.fq") || die "cant make temporary fastq file";

# read in manta records
while(<VCF>){
    chomp;
    next if /^#/;
    my @F = split("\t",$_);

    # get the contig tag, if it exists, and print sequence to a temp fastq file
    if ($F[7] =~ /CONTIG=([ACGTNactgn]+)/){
	my $contig = $1;
	my $quals = '#' x length($contig);
	print FQ "\@$F[2]\n$contig\n+\n$quals\n";
    }
}
close VCF;
close FQ;

my %hits = (); # hash of hits
my %totalhits = ();

# map fastq file to reference with minimap2 and iterate through records
open(SAM,"minimap2 -N 50 -p 0.5 --mask-level 0.8 -ax sr $refseq $tmp/input.fq | $tee $tmp/minimap.sam |") || die;
open(D,">$tmp/minimap.tsv") || die;
while(<SAM>){
    chomp;

    next if /^@/;
    
    my @F = split("\t",$_);

    $totalhits{$F[0]}++;
      
    my $flag = $F[1];

    next if $flag & hex("0x800");  #skip if its a supp alignment
    
    next if !/SA:Z/;  # must have a supplementary alignment in the SA tag

    # get supplemental hit
    /SA:Z:(\S+);/;
    my @l=split(",",$1); 

    # check strand of primary hit
    my $pstrand = "+"; 
    $pstrand = "-" if $flag & hex("0x10"); 

    /NM:i:(\d+)/;
    my $pnm = $1;
    
    # get the alignment coordinates on the read from the cigar string (e.g., positions X through Y are aligned)
    my $cigar = $F[5];
    my $st = 0; # align start
    my $alen = 0; # align length
    my $qlen = 0; # read length
    while($cigar =~ /(\d+)([MSHDI])/g){ # get cigar operation
	my $s = $1;
	my $o = $2;
	$qlen += $s if $s ne 'D'; # add to the query length unless its a del
	$st += $s if ($o =~ /S|H/ and $alen == 0); # change start position iff this is the first operation and its a clip 
	$alen += $s if ($o =~ /[MI]/); # add to the alignment length if its a M or I
    }
    my $pid = 1 - $pnm / $alen;
    my @p = ($st+1,$st+$alen-1); # make array with (start,end)
    @p = reverse(map { $qlen - $_ } @p) if ($pstrand eq '-'); # reverse if its a minus strand hit

    map { die "parse cigar failed" if $_ < 0 or $_ > $qlen } @p;
    
    # do same for supplemental cigar
    my $scigar = $l[3];
    my $snm = $l[5];
    my $sstrand = $l[2];
    $st = 0;
    $alen = 0;
    $qlen = 0;
    while($scigar =~ /(\d+)([MSHDI])/g){
        my $s = $1;
        my $o = $2;
        $qlen += $s if $s ne 'D'; # add to the query length unless its a del
        $st += $s if ($o =~ /S|H/ and $alen == 0); # change start position iff this is the first operation and its a clip
	$alen += $s if ($o =~ /[MI]/); # add to the alignment length if its a M or I
    }
    my $sid = 1 - $snm / $alen;
    my @s = ($st,$st+$alen-1);
    @s = reverse(map { $qlen - $_ } @s) if ($sstrand eq '-');

    map { die "parse cigar failed" if $_ < 0 or $_ > $qlen } @s;
    
    # calculate MASKLEVEL (1.0 - MASKLEVEL) == fraction of query that overlaps between primary and supplemental hits 
    my $ov = 0;
    unless ($s[0] > $p[1] or $p[0] > $s[1]){
	$ov = (min($p[1],$s[1]) - max($p[0],$s[0])) / min($p[1]-$p[0],$s[1]-$s[0]);
    }

    die "masklevel failed" if $ov < 0 or $ov > 1;
    
    my $mq = ($F[4] < $l[4] ? $F[4] : $l[4]);
    # skip if masklevel threshold exceeded
    next if $ov > (1 - $masklevel);
    
    next unless $mq >= $minMQ && $pid > $fracIdent && $sid > $fracIdent;

    # store hits
    push @{$hits{$F[0]}}, [ $F[2],$F[3]-1,$F[3],$l[0],$l[1]-1,$l[1],$F[0],$mq,$pstrand,$sstrand,$F[5] . ";" . $l[3] ];
    print D join("\t",$F[2],$F[3]-1,$F[3],$l[0],$l[1]-1,$l[1],$F[0],$mq,$pstrand,$sstrand,$F[5] . ";" . $l[3]),"\n";
}
close SAM;
close D;

die "no minimap hits!" if scalar keys %hits == 0;

# get hotspot annotations using bedtools. store in hash by mantaID
my %knowntrans = ();
open(BT,"$svtools vcftobedpe -i $ARGV[0] 2> /dev/null | $bedtools pairtopair -is -slop $slop -a stdin -b $knownfile |") || die "cant run bedtools";
while(<BT>){
    chomp;
    my @F = split("\t",$_);

    # get orientation of manta hit
    my $orientation= 'same';
    $orientation = 'opposite' if ($F[14] =~ /^[ACTGactg]+\]|\[[ACTGactg]+$/);

    # store name if right orientation
    if ($F[$#F] eq '.' or ($orientation eq 'same' && $F[$#F] eq $F[$#F-1]) or ($orientation eq 'opposite' && $F[$#F] ne $F[$#F-1])){
	$knowntrans{$F[12]} = $F[$#F-3];
	$knowntrans{$F[15]} = $F[$#F-3];
    }
}
close BT;

# outfile
open(O,">$out") || die;

#open manta VCF again
open(VCF,"zcat $ARGV[0] |") || die "Cant open manta vcf file $ARGV[0]";
while(<VCF>){
    # add info tag and filters
    if (/^##FILTER/){
	print O '##INFO=<ID=KNOWNSV,Number=.,Type=String,Description="Known hotspot translocation">',"\n";
	print O '##INFO=<ID=CONTIGHITS,Number=.,Type=Integer,Description="Number of high-quality primary genomic alignments to manta contig">',"\n";
	do {
	    print O;
	    $_ = <M>;
	} while(/^##FILTER/);
	print O '##FILTER=<ID=LowReads,Description="Failed minimum number of PR or SR reads">',"\n";
	print O '##FILTER=<ID=FailedContig,Description="Mapping of contig failed to reproduce breakends">',"\n";
	print O '##FILTER=<ID=FailedCov,Description="Coverage of DEL/DUP does not support an SV call">',"\n";
	next;
    } elsif (/^#/){
	print O;
	next;
    }
    
    # handle records
    chomp;
    my @F = split("\t",$_);

    my @fmtk = split(":",$F[8]);
    my @fmtv = split(":",$F[9]);

    my %fmt = ();
    map { $fmt{$fmtk[$_]} = $fmtv[$_]; $fmt{$fmtk[$_]} = [ split(",",$fmt{$fmtk[$_]}) ] if $fmt{$fmtk[$_]}=~/,/; } 0..$#fmtk;

    $F[7] =~ /SVTYPE=([^;]+)/;
    my $svtype = $1;
    
    # add hotspot tag, if exists
    $F[7] .= ";KNOWNSV=" . $knowntrans{$F[2]} if (defined($knowntrans{$F[2]}));	           
    
    # add low reads filter if PR and SR reads are low or VAF is too low
    my $svabund = 0.0;
    my @filter = ();
    @filter = split(',',$F[6]) if $F[6] ne 'PASS';
    
    if (defined($fmt{PR}) and defined($fmt{SR}) and (max(@{$fmt{PR}}) + max(@{$fmt{SR}})) > 0){
	$svabund = (($fmt{PR}->[1] + $fmt{SR}->[1]) / ($fmt{PR}->[0] + $fmt{PR}->[1] + $fmt{SR}->[0] + $fmt{SR}->[1])) * 100;
    }
    
    if (!defined($fmt{PR}) or !defined($fmt{SR}) or $fmt{PR}->[1] < $PR or $fmt{SR}->[1] < $SR or $svabund < $minSVabund){
	push @filter, "LowReads";
	
    } elsif (defined($fmt{DHFFC}) and defined($fmt{DHBFC}) and ($svtype =~ /DEL/ and $fmt{DHFFC} > $delcovratio) or ($svtype =~ /DUP/ and $fmt{DHBFC} < $dupcovratio)){
	push @filter, "FailedCov";
	
	# add no contig filter if no contig/imprecise breakends
    } elsif ($F[7] !~ /CONTIG=/ or !defined($hits{$F[2]})){
	push @filter, "FailedContig";
	
    } else {
	
	my $orientation= 'same';
	my $chr1 = '';
	my $chr2 = '';
	my $pos1 = '';
	my $pos2 = '';
	
	# get positions
	if ($F[2] =~ /DEL|INV|DUP|INS/){
	    $F[7] =~/END=(\d+)/;
	    
	    $chr1 = $F[0];
	    $pos1 = $F[1];
	    $chr2 = $F[0];
	    $pos2 = $1;
	    
	} elsif ($F[2] =~ /BND/) {
	    $F[4] =~/([^\[\]]+):(\d+)/;
	    
	    $chr1 = $F[0];
	    $pos1 = $F[1];
	    $chr2 = $1;
	    $pos2 = $2;
	    
	    $orientation = 'opposite' if ($F[4] =~ /^[ACTGactg]+\]|\[[ACTGactg]+$/);	    
	}

	# add hotspot tag, if exists
	$F[7] .= ";CONTIGHITS=" . $totalhits{$F[2]};
	
	# get hits that overlap, allowing for some $slop
	my $foundhit = 0;
	foreach my $i (@{$hits{$F[2]}}){
	    $foundhit = 1 if ((($i->[0] eq $chr1 && $pos1 < $i->[2]+$slop && $pos1 > $i->[1]-$slop &&
			      $i->[3] eq $chr2 && $pos2 < $i->[5]+$slop && $pos2 > $i->[4]-$slop) ||
			      ($i->[0] eq $chr2 && $pos2 < $i->[2]+$slop && $pos2 > $i->[1]-$slop &&
			       $i->[3] eq $chr1 && $pos1 < $i->[5]+$slop && $pos1 > $i->[4]-$slop)) &&
			       (($orientation eq 'same' && $i->[8] eq $i->[9]) || ($orientation eq 'opposite' && $i->[8] ne $i->[9])));
		}
	push @filter, "FailedContig" if $foundhit == 0;
	
    }
    $F[6] = join(';',@filter);
    $F[6] = "PASS" if scalar @filter == 0;
    
    print O join("\t",@F),"\n";
}
close VCF;
close O;    

exit;
