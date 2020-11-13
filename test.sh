#!/bin/bash

output=test-output.genozip

is_windows=`uname|grep -i mingw`
is_mac=`uname|grep -i Darwin`

hg19=data/hg19.p13.plusMT.full_analysis_set.ref.genozip
GRCh38=data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip

arg1=$1

# debug
is_debug=`echo $1|grep debug`
if [ -n "$is_debug" ]; then 
    debug=-debug; 
    shift
fi

# -----------------
# platform settings
# -----------------
if [ -n "$is_windows" ]; then
    genozip=./genozip${debug}.exe
    genounzip=./genounzip${debug}.exe
    genocat=./genocat${debug}.exe
    genols=./genols${debug}.exe
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    genozip=./genozip${debug}
    genounzip=./genounzip${debug}
    genocat=./genocat${debug}
    genols=./genols${debug}
    path=$PWD/
fi

exes=($genozip $genounzip $genocat $genols)
for exe in ${exes[@]}; do
    if [ ! -x $exe ]; then
        echo "Error: $exe does not exist"
        exit 1
    fi
done

if `command -v md5 >& /dev/null`; then
    md5=md5 # mac
else
    md5=md5sum 
fi

cmp_2_files() {
    if (( `$md5 $1 ${2%.*} | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then
        echo "MD5 comparison FAILED:"
        $md5 $1 ${2%.*}
        exit 1
    fi
}

test_header() {
    sep="=======================================================================================================\n"
    printf "\n${sep}TESTING $1 \n${sep}"
}

test_count_genocat_lines() {
    local cmd="$genocat $output $2"
    test_header "$cmd"
    $genozip $arg1 -fo $output || exit 1
    local wc=`$cmd | wc -l`
    if (( $wc != $3 )); then
        echo "FAILED - expected $3 lines, but getting $wc"
        exit 1
    fi  
}

test_bam() {
    if `command -v samtools >& /dev/null`; then
        test_header "$1 - input and output as BAM"
        grep -v TooBigForSamTools $1 > bam-test.input.sam || exit 1
        local arg=("--no-PG" "")
        samtools view --help |& grep no-PG >& /dev/null # test for --no-PG (exists since samtools 1.10, see https://github.com/samtools/samtools/releases/)
        local missing_no_PG=$? # $?==0 if exists, 1 if not
        samtools view ${arg[$missing_no_PG]} bam-test.input.sam -OBAM -h > bam-test.input.bam || exit 1
        $genozip $arg1 bam-test.input.bam $2 -fto $output || exit 1
        $genounzip $arg1 $output $2 --force --output bam-test.output.bam || exit 1
        sleep 0.2 # wait for BAM to be flushed to the disk
        cmp_2_files bam-test.input.bam bam-test.output.bam.fake-extension-removed-by-cmp_2_files
        rm -f bam-test.input.sam bam-test.input.bam bam-test.output.bam
    fi
}

# minimal files - expose edge cases where fields have only 1 instance
files=(minimal.vcf minimal.sam minimal.fq minimal.fa minimal.gvf minimal.genome_Full.me23.txt)
for file in ${files[@]}; do
    test_header "$file - minimal file test"

    if [ ! -f $file ] ; then echo "$file: File not found"; exit 1; fi

    cat $file | tr -d "\r" > unix-nl.$file || exit 1
    $genozip $arg1 unix-nl.$file -ft -o $output || exit 1
    rm -f unix-nl.$file $output
done

test_bam test-file.sam

files=(test-file.vcf test-file.sam test-file.fq test-file.fa test-file.gvf test-file.genome_Full.me23.txt)
for file in ${files[@]}; do
    test_header "$file - basic test - Unix-style end-of-line"

    if [ ! -f $file ] ; then echo "$file: File not found"; exit 1; fi

    cat $file | tr -d "\r" > unix-nl.$file
    $genozip unix-nl.$file -ft -o $output || exit 1

    test_header "$file - Window-style end-of-line"

    if [ ! -n "$is_mac" ]; then  # note: sed on mac doesn't recognize \r
        sed 's/$/\r/g' unix-nl.$file > windows-nl.$file || exit 1 # note: sed on mac doesn't recognize \r
        $genozip $arg1 windows-nl.$file -ft -o $output || exit 1
        rm -f unix-nl.$file windows-nl.$file
    fi

    test_header "$file - as URL"
    $genozip $arg1 file://${path}$file -ft -o $output || exit 1

    test_header "$file - encrypted"
    $genozip $arg1 $file --password abc -ft -o $output || exit 1

    test_header "$file - redirected from stdin"
    cat $file | $genozip $arg1 --test --force --output $output --input-type ${file#*.} - || exit 1

    if [ $file != test-file.sam ] && [ $file != test-file.genome_Full.me23.txt ]; then
        allow_compressed=1;
    else 
        allow_compressed=0;
    fi

    if `command -v gzip >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with gzip"
        cp -f $file copy.$file
        gzip -f copy.$file
        $genozip $arg1 copy.${file}.gz -ft -o $output || exit 1
        rm -f copy.${file}.gz
    fi
    
    if `command -v bzip2 >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with bzip2"
        
        # bzip2 fails if it can't chmod the dst file, but chmod from WSL fails on an NTFS filesystem
        if [ -n "`uname -a | grep -i microsoft-standard`" ]; then
            copy=~/copy.$file
        else 
            copy=./copy.$file
        fi

        cp -f $file $copy
        bzip2 $copy
        $genozip $arg1 ${copy}.bz2 -ft -o $output || exit 1
        rm -f ${copy}.bz2
    fi
    
    if `command -v xz >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with xz"
        cp -f $file copy.$file
        xz -f copy.$file
        $genozip $arg1 copy.${file}.xz -ft -o $output || exit 1
        rm -f copy.${file}.xz
    fi
        
    if [ -z "$is_windows" ]; then # windows can't redirect binary data
        test_header "$file - redirecting stdout"
        $genozip $arg1 ${file} --stdout > $output || exit 1
        $genounzip $arg1 $output -f || exit 1
        cmp_2_files $file $output
    fi

    test_header "$file - non-bound multiple files"
    file1=copy1.$file
    file2=copy2.$file
    cp -f $file $file1
    cp -f $file $file2
    $genozip $arg1 $file1 $file2 -ft || exit 1
    $genounzip $arg1 ${file1}.genozip ${file2}.genozip -t || exit 1
    rm -f $file1 $file2 ${file1}.genozip ${file2}.genozip

    test_header "$file - bind & unbind"
    file1=copy1.$file
    file2=copy2.$file
    cp -f $file $file1
    cat $file | sed 's/PRFX/FILE2/g' > $file2
    $genozip $arg1 $file1 $file2 -ft -o $output || exit 1
    $genounzip $arg1 $output -u -t || exit 1
    rm -f $file1 $file2

    test_header "$file --optimize - NOT checking correctness, just that it doesn't crash"
    $genozip $arg1 $file -f --optimize -o $output || exit 1

done

# files represent cases that cannot be fit into test-file.* because they would conflict
files=(test-file-domqual.fq test-file-domqual.sam test-file-unaligned.sam)
for file in ${files[@]}; do

    if [ ! -f $file ] ; then echo "$file: File not found"; exit 1; fi

    test_header "$file - special case test"
    cat $file | tr -d "\r" > unix-nl.$file
    $genozip $arg1 unix-nl.$file -ft -o $output || exit 1
done

# Test SAM reconstruction of BAM files
#files=(td/test.NA12878.chr22.1x.bam td/test.m64136_200621_234916.ccs.10k.bam)
for file in ${files[@]}; do
    if [ ! -f $file ] ; then echo "$file: File not found"; exit 1; fi

    test_header "$file - genocat bam.genozip and compare to original SAM"

    $genozip $arg1 -f bug.bam -o $output
    $genocat $output > copy.sam 
    cmp_2_files $file copy.sam
    rm -f $copy.sam $output
done

# Test binding SAM files with lots of contigs (no reference)
test_header "binding SAM files with lots of contigs (no reference)"
file=td/test.transfly-unsorted.sam
cp -f $file copy.unsorted1.sam
cp -f $file copy.unsorted2.sam
$genozip $arg1 copy.unsorted1.sam copy.unsorted2.sam -ft -o $output || exit 1
rm -f copy.unsorted1.sam copy.unsorted2.sam $output

# VCF gtshark test
if `command -v gtshark >& /dev/null`; then
    test_header "test-file.vcf --gtshark"
    $genozip $arg1 test-file.vcf --gtshark -ft -o $output || exit 1
fi

# FASTA genocat tests
test_count_genocat_lines test-file.fa "--sequential" 9
test_count_genocat_lines test-file.fa "--header-only" 3
test_count_genocat_lines test-file.fa "--header-one" 3
test_count_genocat_lines test-file.fa "--no-header" 15
test_count_genocat_lines test-file.fa "--no-header --sequential" 6
test_count_genocat_lines test-file.fa "--grep cytochrome" 6
test_count_genocat_lines test-file.fa "--grep cytochrome --sequential " 2
test_count_genocat_lines test-file.fa "--grep cytochrome --sequential --no-header " 1

# FASTQ genocat tests
test_count_genocat_lines test-file.fq "--header-only" `grep @ test-file.fq | wc -l` 
test_count_genocat_lines test-file.fq "--header-one" `grep @ test-file.fq | wc -l`
test_count_genocat_lines test-file.fq "--grep line5" 4
test_count_genocat_lines test-file.fq "--grep line5 --header-only" 1

#files=`ls backward-compatibility-test/*.genozip` 
#for file in $files; do
#    test_header "$file - backward compatability test"
#
#    if [ `basename $file .vcf.genozip` = test-file.1.1.3 ]; then # in v1 we didn't have the -t option
#        $genounzip $arg1 ${file} -fo $output || exit 1
#        cmp_2_files backward-compatibility-test/test-file.1.1.3.vcf $output
#    else
#        $genounzip $arg1 -t $file || exit 1
#    fi
#done

test_header "test-file.vcf without FORMAT or samples"
file=test-file.vcf
cut -f1-8 $file > copy.$file
$genozip $arg1 copy.$file -ft -o $output || exit 1
rm -f copy.$file

test_header "subsets (~3 VBs) or real world files"
rm -f td/*.genozip
$genozip $arg1 -ft td/*.sam* td/*.vcf* td/*.fa* td/*.f*q* td/*gvf* td/*genome* || exit 1

test_header "--make-reference"
file=test-file-ref.fa 
$genozip $arg1 --make-reference $file --force -o copy.${file}.ref.genozip || exit 1

test_header "unaligned SAM with --reference"
$genozip $arg1 -f --md5 --reference copy.${file}.ref.genozip test-file-unaligned.sam --test || exit 1

test_header "unaligned SAM with --reference - from stdin"
cat test-file-unaligned.sam | $genozip $arg1 -f --md5 --reference copy.${file}.ref.genozip --test --input sam --output $output - || exit 1

test_bam test-file-unaligned.sam -ecopy.${file}.ref.genozip

rm -f copy.${file}.ref.genozip $output

test_header "command line with mixed SAM and FASTQ files with --reference"
echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
rm -f td/*.genozip
$genozip $arg1 -f --md5 --reference $GRCh38 td/test.transfly-unsorted.sam td/test.transfly.fq td/test.transfly-sorted.sam || exit 1
$genounzip $arg1 -t -e $GRCh38 td/test.transfly-unsorted.sam.genozip td/test.transfly-sorted.sam.genozip || exit 1

test_header "multiple bound SAM with --REFERENCE" 
rm -f td/*.genozip
$genozip $arg1 -f --md5 --REFERENCE $GRCh38 td/test.transfly-unsorted.sam td/test.transfly-sorted.sam -o $output || exit 1
$genounzip $arg1 -t $output || exit 1

test_header "SAM with --reference and --password" 
rm -f td/*.genozip
$genozip $arg1 -f --md5 --reference $GRCh38 td/test.transfly-unsorted.sam --password 123 -o $output || exit 1
$genounzip $arg1 -t --reference $GRCh38 -p 123 $output || exit 1

test_header "SAM with --REFERENCE and --password" 
rm -f td/*.genozip
$genozip $arg1 -f --md5 --REFERENCE $GRCh38 td/test.transfly-unsorted.sam --password 123 -o $output || exit 1
$genounzip $arg1 -t -p 123 $output || exit 1

test_header "paired FASTQ with --reference, --password and --md5"
rm -f td/*.genozip
$genozip $arg1 -f --md5 -e $GRCh38 --pair -p 1234 td/test.divon-R1.100K.fq.bz2  td/test.divon-R2.100K.fq.bz2 -o td/pair.genozip || exit 1
$genounzip $arg1 -t --password 1234 -e $GRCh38 td/pair.genozip || exit 1

test_header "4 paired FASTQ with --REFERENCE"
rm -f td/*.genozip
file1=td/test.divon-R1.100K.fq.bz2
file2=td/test.divon-R2.100K.fq.bz2
file3=copy.$(basename $file1)
file4=copy.$(basename $file2)
cp -f $file1 $file3
cp -f $file2 $file4
$genozip $arg1 --force -m2E $GRCh38 --output $output $file1 $file2 $file3 $file4 || exit 1
$genounzip $arg1 -t $output || exit 1
rm -f $file3 $file4 $output

test_header "multiple bound VCF with --reference (hg19), and unbind"
rm -f td/*.genozip
file1=copy1.test-file.vcf
file2=copy2.test-file.vcf
cp -f td/test.GFX0241869.filtered.snp.vcf $file1
cp -f td/test.GFX0241869.filtered.snp.vcf $file2
$genozip $arg1 -f --md5 --reference $hg19 $file1 $file2 --output $output || exit 1
$genounzip $arg1 -t -e $hg19 --unbind $output || exit 1
rm -f $file1 $file2 $output

test_header "multiple VCF with --REFERENCE using hg19" 
rm -f td/*.genozip
$genozip $arg1 -f --md5 --REFERENCE $hg19 td/test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf td/test.GFX0241869.filtered.snp.vcf || exit 1
$genounzip $arg1 -t td/test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.genozip td/test.GFX0241869.filtered.snp.vcf.genozip || exit 1

printf "\nALL GOOD!\n"

rm -f $output $output td/*.genozip
