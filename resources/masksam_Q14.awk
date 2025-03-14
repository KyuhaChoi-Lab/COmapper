#!/usr/bin/awk -f

################################
#" #!: shebang; specifies an interpreter for the instructions in a script
#" /usr/bin/awk: the interpreter
#" -f: to specify a file that contains awk script
################################
# filter-out not-aligned reads before putting as an input

# sam field
# 1: read id
# 2: flag
# 3: ref
# 4: pos
# 5: mapq
# 6: cigar
# 7: rnext
# 8: pnext
# 9: tlen
# 10: seq
# 11: qual
#
# Add header from the original bam file after masking
#
# P-error cutoff: Q=14 ; P_error=0.03981

BEGIN{FS="\t"; OFS="\t"}

/^@/{next}

{ if ( $5 == 0 ) next;}
{ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, $3, $4, $5, $6, $7, $8, $9 }
{
    n = split($10, seq, "")
    n = split($11, qual, "")
    for (i = 1; i <= n; i++) {
        # ASCII 33-47 corresponds to Q=0 to Q=14
        if (qual[i] ~ /[!"#$%&'()*+,-./]/)
            printf "N"
        else
            printf "%s", seq[i]
    }
    printf "\t%s\n", $11
}
