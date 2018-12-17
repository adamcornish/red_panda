#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Basic qw(:all);
use Getopt::Std;

sub rm_variants_with_low_read_count;
sub rm_homozygous_variants;
sub create_directories;
sub printtime;
sub print_debug_statement;

$"="\n";

my %opt;
getopts ( "hq:v:s:d:g:A:V", \%opt );
die "Usage: $0 -q quant.sf -v variants.vcf -s sample_name -d output_directory -g (hg38|mm10) -A /path/to/genome/annotation.gtf -V" if $opt{h};

print "genome: $opt{g}\n";
die "You must provide a sample name.\n" unless $opt{s};
die "You must provide a directory for output.\n" unless $opt{d};
die "Directory $opt{d} does not exist.\n" unless -d $opt{d};
die "Only genomes hg38 and mm10 are supported.\n" if ( $opt{g} ne "hg38" and $opt{g} ne "mm10" );
die "No annotation file provided.\n" unless $opt{A};

# create directories if they don't exist
create_directories( $opt{d} );

# Take the list of isoforms previously identified as expressed and create a bed file consisting of only those isoforms
#     Generate a bed file for each isoform
#     Intersect the VCF file with the bed file
# Only grab isoforms with a TPM >= 1 (only works with mm10 and hg38 (probably other human ones, but untested))
#     Search for the ENST ID to get the ENSG ID. 
#     Then search for the ENSG ID to get ALL of the ENST IDs for that gene. 
#     Then loop through THAT list to see if any of them are in the expressed quant list.

# TODO turn this into a function
if ( $opt{q} )
{
    chomp ( my @salmon_transcripts = `grep -vP "\\S+\\s\\S+\\s\\S+\\s0" $opt{q} | grep -v Name | sort` );
    my $annotation_file = $opt{A}; # "data/annotation/ref-transcripts.gtf"; # grab this from bcbio
    my $do_not_check    = "";
    my $ENST_string     = "ENST";
    my @ENST_isoforms;
    $ENST_string = "ENSMUST" if $opt{g} eq "mm10";

    # build hash with transcript ids as keys
    printtime();
    print "Build hashes with transcript IDs as keys\n";

    my %ENST_hash;
    my %ENST_hash_full_entry;
    open IN, "<$annotation_file";
    while ( my $line = <IN> )
    {
        my ($key) = $line =~ /($ENST_string.+?)"/;
        if ( $line =~ /transcript_biotype "protein_coding/ )
        {
            $ENST_hash{$key} = $line if $line =~ /\stranscript\s/;
            $ENST_hash_full_entry{$key} .= $line;
        }
    }
    close IN;

    # build hash with gene ids as keys
    printtime();
    print "Build hash with gene IDs as keys\n";
    my %ENSG_hash;
    open IN, "<$annotation_file";
    while ( my $line = <IN> )
    {
        if ( $line =~ /\stranscript\s/ )
        {
            my ($key) = $line =~ /gene_id "(.+?)"/;
            if ( $line =~ /transcript_biotype "protein_coding"/ )
            {
                $ENSG_hash{$key} = $ENSG_hash{$key} ? "$ENSG_hash{$key}!$line" : $line;
            }
        }
    }
    close IN;
    printtime();
    print "Finished building hashes\n";

    # Get a list of exons to search through by pushing lines from the annotation file into @ENST_isoforms
    my $counter = 0;

    printtime();
    print "Generate an annotation file containing only transcripts that are expressed (TPM > 1) using $opt{q} as input.\n";
    foreach my $line ( @salmon_transcripts )
    {
        if ( $opt{V} )
        {
            printtime();
            printf ( "Working (%04d / %04d): %s\n", $counter++, $#salmon_transcripts+1, $line );
        }

        my ($ENST)     = $line =~ /($ENST_string.+?)\s/;
        unless ( $do_not_check =~ /$ENST/ )
        {
            my $tmp = $ENST =~ /$do_not_check/;

           #my $ENSG_line  = `grep $ENST $annotation_file | grep -P "\\stranscript\\s"`; # old line
            my $ENSG_line = $ENST_hash{$ENST};

            if ( $ENSG_line )
            {
                my ($ENSG)     = $ENSG_line =~ /gene_id "(.+?)"/;
                 
                # Get all of the transcripts that carry this gene id AND are protein coding transcripts
               #my @ENST_lines = `grep $ENSG $annotation_file | grep -P "\\stranscript\\s" | grep 'transcript_biotype "protein_coding'`; # old line
                my @ENST_lines = split /!/, $ENSG_hash{$ENSG};

                if ( $#ENST_lines > -1 ) # make sure there are actually protein coding transcripts for this gene
                {
                    my ($isoforms) = $ENSG_line =~ /gene_name "(.+?)"/;
                    my $DNC_list   = ""; #DNC = Do Not Check

                    $isoforms .= ",$ENSG:";

                    foreach my $E_line ( @ENST_lines )
                    {
                        my ($E)    = $E_line =~ /($ENST_string.+?)"/;
                        my $exists = grep ( /$E/, @salmon_transcripts );
                        $isoforms .= "$E;" if $exists;
                        $DNC_list .= "$E;" if $exists;
                    }
                    $do_not_check .= "$DNC_list:";
                    chop $DNC_list;
                    push @ENST_isoforms, $isoforms;
                }
            }
        }
    }

    # @ENST_isoforms will have a list of all the isoforms currently expressed, so lets dump all the records out to a gff file
    my ($sample_name) = $opt{s};
    open OUT, ">$opt{d}/quant/$sample_name.gff";
    printtime();
    print "Printing expressed isoforms out to $opt{d}/quant/$sample_name.gff\n";
    foreach my $line ( @ENST_isoforms )
    {
        my @isoforms = $line =~ /($ENST_string\d+)/g;
        foreach my $isoform ( @isoforms )
        {
            # only print out the exon lines, because including the full entry means that variants can be called in the introns
            my @lines = split /\n/, $ENST_hash_full_entry{$isoform};
            my @exon_lines = grep /\sexon\s/, @lines;
            print OUT "@exon_lines\n";
        }
    }
    close OUT;
}



# TODO turn this part into a function
if ( $opt{v} )
{
    my ($sample_name) = $opt{s};
    open GOOD_VAR,     ">$opt{d}/vcf/good/$sample_name.vcf";
    open BAD_VAR,      ">$opt{d}/vcf/bad/$sample_name.vcf";
    open TO_CHECK_VAR, ">$opt{d}/vcf/to_check/$sample_name.vcf";
    my $head = `grep "^#" $opt{v}`;
    print GOOD_VAR $head;
    print BAD_VAR $head;
    print TO_CHECK_VAR $head;
    my @good_variants;
    my @bad_variants;
    my @to_check_variants;
    chomp ( my @unfiltered_variants = `cat $opt{v} | grep -v "^#"` ); 
    printtime();
    print "Total number of unfiltered : $#unfiltered_variants\n";

    # remove anything that only has one read supporting it
    my @variants;
    foreach my $var ( @unfiltered_variants )
    {
        my ($DP,$ref,$alts) = $var  =~ /DP=(\d+).+:(\d+),(.+)/;
        my @alt         = $alts =~ /(\d+)/g;
        my $alt_sum     = 0;
        foreach my $item ( @alt ) { $alt_sum += $item; }
        push @variants, $var if $alt_sum > 1 and $DP > 9;
    }

    printtime();
    print "Total variants w/o alt=1   : $#variants\n";

    my @pct_filtered_variants = &rm_variants_with_low_read_count ( @variants );
    printtime();
    print "Filtering out variants with low read counts : $#pct_filtered_variants\n";

    # this array contains variants that AREN'T homozygous-looking
    my @hom_filtered_variants = &rm_homozygous_variants ( $sample_name, @pct_filtered_variants );
    printtime();
    print "Filtering out homozygous-looking variants : $#hom_filtered_variants\n";
    
    open OUT, ">>$opt{d}/vcf/filtered/$sample_name.vcf";
    system "grep -P '^#' $opt{v} > $opt{d}/vcf/filtered/$sample_name.vcf";
    print OUT "@hom_filtered_variants";
    close OUT;

    # Get the intersection of the gff file containing the isoforms that are expressed and the vcf file that has been filtered up to this point.
    # Remove all of the isoforms from our previous list that do not have putative heterozygous variants in them.
    printtime();
    print "Start intersect of $sample_name based on expression\n";

    # Grab the header of the vcf file and only keep variants that are in the boundaries of exons
    system "grep \"^#\" $opt{v} > $opt{d}/vcf/expressed/$sample_name.vcf";
    system "bedtools intersect -a $opt{d}/vcf/filtered/$sample_name.vcf -b $opt{d}/quant/$sample_name.gff -u >> $opt{d}/vcf/expressed/$sample_name.vcf";

    printtime();
    print "End intersect of $sample_name based on expression\n";

    # Build hash containing each isoform that was expressed where $key = transcript ID (e.g., ENST[\d+] or ENMUST[\d+])
    chomp ( my @isoforms = `cat $opt{d}/quant/$sample_name.gff` );
    my %tsids;

    foreach my $line ( @isoforms )
    {
        my ($gsym) = $line =~ /transcript_id "(.+?)"/;
        $tsids{$gsym} = "" unless $tsids{$gsym};
        $tsids{$gsym} .= "$line\n";
    }

    my $hash_size = keys %tsids; # total isoforms to search through
    my $i = 1;                   # counter to track progress

    # loop through each isoform and get the variants present in them and perform the bimodal expression tests
    foreach my $key ( keys %tsids )
    {
        # TODO this is janky as heck; fix it so we aren't writing to disk for every transcript being checked
        open OUT, ">.$sample_name.tmp.gff";
        print OUT $tsids{$key};
        close OUT;
        my ($gene_id) = $tsids{$key} =~ /gene_name "(.+?)"/;

        if ( $opt{V} )
        {
            printtime();
            printf ( "Working on %s (%04d / %04d): %s\n", $key, $i++, $hash_size, $gene_id );
        }

        chomp ( my @vars = `bedtools intersect -a $opt{d}/vcf/expressed/$sample_name.vcf -b .$sample_name.tmp.gff -u` );
        my @high_fractions;
        my @low_fractions;
        my $max_fraction = 0;

        if ( $#vars > -1 )
        {
            foreach my $var ( @vars )
            {
                my ($pos,$DP,$ref,$alt) = $var  =~ /(chr\S+\s\d+).+DP=(\d+).+:(\d+),(.+?)(?:,|$)/;
                my $fraction = $alt/$DP;

                # Create two bins: 
                #     a.Alt allele fraction > 50% (major)
                #     b.Alt allele fraction < 50% (minor)
                push @high_fractions, "$key|$gene_id|$pos|$DP|$fraction|$var" if $fraction >= 0.5;
                push @low_fractions,  "$key|$gene_id|$pos|$DP|$fraction|$var" if $fraction < 0.5;
                $max_fraction = $fraction if $fraction > $max_fraction; # keep track of the highest fraction calculated
            }

            my $bimodal = 0; # is the distribution bimodal?
            my $bia_hfrac = 0;
            my $bia_lfrac = 0;

            # Search through the bins to determine if the alt alleles agree with a bimodal distributation, and if they do, simply add these variants to the good_variants list
            # TODO: include some logic where the depth of the variant in question adds weight to the decision-making process
            foreach my $high_item ( @high_fractions )
            {
                my ($hkey,$hgene_id,$hpos,$hDP,$hfrac,$hvcf_line) = split /\|/, $high_item;

                # If we already have a bimodal distribution established and this fraction is within an acceptable range, add it to the good_variants list
                if ( $hfrac > ( $bia_hfrac - 0.05 ) and $hfrac < ( $bia_hfrac + 0.05 ) )
                {
                    push @good_variants, $hvcf_line;
                }
                else
                {
                    my $found = 0; # did we find a bimodal distribution?

                    # Start with the list of major alleles and determine if its fractions are tightly clustered;
                    #   if so, any valid minor alleles should be 1-major allele value. Start by evaluating the alt alleles with the highest fraction of reads supporting it.
                    for ( my $i = 0; $i < @low_fractions; $i++ )
                    {
                        my $low_item = $low_fractions[$i];
                        my ($lkey,$lgene_id,$lpos,$lDP,$lfrac,$lvcf_line) = split /\|/, $low_item;
                        my $sum = $hfrac + $lfrac;
                        if ( $sum > 0.90 and $sum < 1.1 )
                        {
                            push @good_variants, $hvcf_line if !$found; # only add the hvcf_line once
                            $bimodal = 1;
                            $found     = 1;
                            $bia_hfrac = $hfrac;
                            $bia_lfrac = $lfrac;
                            push @good_variants, $lvcf_line;
                        }
                        else
                        {
                            # If there exist variants that do not agree with the bimodal distribution, place these in a false_positive list.
                            push @bad_variants, $lvcf_line;
                        }
                    }

                    # If there are not a sufficient number of variants present in this gene or there is no obvious way to determine a bimodal distribution,
                    # place these variants in the list of variants to be evaluated using the consensus method
                    push @to_check_variants, $hvcf_line if !$found;
                }
            }

            foreach my $low_item ( @low_fractions ) 
            {
                my ($lkey,$lgene_id,$lpos,$lDP,$lfrac,$lvcf_line) = split /\|/, $low_item;
                
                # If there are not a sufficient number of variants present in this gene or there is no obvious way to determine a bimodal distribution,
                # please these variants in the list of variants to be evaluated using the consensus method
                $lfrac < 0.1 ? push @bad_variants, $lvcf_line : push @to_check_variants, $lvcf_line;

                # TODO If there are not any major alleles, it’s still possible that the minor alleles present are valid
            }

           #print_debug_statement (\@high_fractions, \@low_fractions, \@good_variants, \@bad_variants, \@to_check_variants );
           
        }
    }

    print GOOD_VAR "@good_variants";
    close GOOD_VAR;
    print BAD_VAR "@bad_variants";
    close BAD_VAR;
    print TO_CHECK_VAR "@to_check_variants";
    close TO_CHECK_VAR;
    system "cp $opt{d}/vcf/good/$sample_name.vcf x && cat x \| vcf-sort -c \| uniq > $opt{d}/vcf/good/$sample_name.vcf";
    system "cp $opt{d}/vcf/bad/$sample_name.vcf x && cat x \| vcf-sort -c \| uniq > $opt{d}/vcf/bad/$sample_name.vcf";
    system "cp $opt{d}/vcf/to_check/$sample_name.vcf x && cat x \| vcf-sort -c \| uniq > $opt{d}/vcf/to_check/$sample_name.vcf";
    system "rm .$sample_name.tmp.gff";
}

# Remove any variants that: have a total read depth < 10, or the fraction of reads supporting the variant does not meet a percentage requirement:
#     a.Fraction must be > 20% if read depthis < 20,
#     b.Fraction must be > 10% if read depth is < 40 and >= 20
#     c.Fraction must be > 5%

sub rm_variants_with_low_read_count
{
    my @variants = @_; 
    my @rval;
    foreach my $var ( @variants )
    {
        my ($ref,$alts) = $var  =~ /.+:(\d+),(.+)/;
        my @alt         = $alts =~ /(\d+)/g;
        my $alt_sum     = 0;
        foreach my $item ( @alt ) { $alt_sum += $item; }
        my $sum         = $alt_sum + $ref;
        if ( $ref ) # don't divide by 0
        {
            my $fraction = $alt_sum / $sum; 

# If the depth is small, use a stricter cutoff for what a minimum fraction of reads is allowed.
#     That is, as we get more reads, then, while the number of reads needed to consider a variant might remain the same (4 reads in this case), the fraction of reads gets smaller
            if    ( $sum <= 20 ) { push @rval, $var if $fraction >= 0.20; }
            elsif ( $sum <= 40 ) { push @rval, $var if $fraction >= 0.10; }
            else                 { push @rval, $var if $fraction >= 0.05; }
        }
        else
        {
            push @rval, $var; # you forgot to put this in initially. you. dumbass.
        }
    }
    return @rval;
}

# TODO add header to the output file
sub rm_homozygous_variants
{
    my $sample_name = shift;
    my @variants = @_;
    my @rval;
    my @hom;
    foreach my $var ( @variants )
    {
        my ($ref,$alt) = $var  =~ /.+:(\d+),(\d+)/;
        my $sum        = $alt + $ref;
        if ( $ref ) # don't divide by 0
        {
            my $fraction = ( $alt / $sum ); 

            # determine if a variant is homozygous-looking based on depth and fraction of reads supporting the variant
            if    ( $sum <= 20  ) { $fraction > 0.86 ? push @hom, $var : push @rval, $var; }
            elsif ( $sum <= 40  ) { $fraction > 0.88 ? push @hom, $var : push @rval, $var; }
            elsif ( $sum <= 60  ) { $fraction > 0.90 ? push @hom, $var : push @rval, $var; }
            elsif ( $sum <= 80  ) { $fraction > 0.92 ? push @hom, $var : push @rval, $var; }
            elsif ( $sum <= 100 ) { $fraction > 0.94 ? push @hom, $var : push @rval, $var; }
            else                  { $fraction > 0.96 ? push @hom, $var : push @rval, $var; }
        }
        else
        {
            push @hom, $var;
        }
    }
    open OUT, ">$opt{d}/vcf/homozygous/$sample_name.vcf";
    print OUT "@hom";
    close OUT;
    return @rval;
}

sub create_directories
{
    my $base_dir = shift;
    system "mkdir -p $base_dir/{vcf,quant,bed}";
    system "mkdir -p $base_dir/vcf/{bad,expressed,filtered,good,homozygous,to_check}";
}

sub printtime
{
    my ($sec,$min,$hour) = localtime(time);
    printf ("[%02d:%02d:%02d] ", $hour, $min, $sec);
}

sub print_debug_statement
{
    # these are all references to the arrays passed to this subroutine and are therefore scalar values
    # @$ dereferences them and allows access to the original array
    my ($high_fractions, $low_fractions, $good_variants, $bad_variants, $to_check_variants) = @_; 
    print "high_fractions   :\n@$high_fractions\n";
    print "low_fractions    :\n@$low_fractions\n";
    print "good_variants    :\n@$good_variants\n";
    print "bad_variants     :\n@$bad_variants\n";
    print "to_check_variants:\n@$to_check_variants\n";
    print "\n", "="x80, "\n\n";
    my $sin = <STDIN>;
}
