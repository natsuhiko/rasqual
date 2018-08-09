#!/bin/bash

if [ ${#@} -le 3 ] ; then
    printf '\nUsage: %s [paired_end|single_end] bam_list_file vcf_input_file vcf_output_file assay_type\n' "${0}" > /dev/stderr;
    printf '\n       assay_type can be \"atac\", \"rna\", or alternatively not specified (general).\n\n' > /dev/stderr;
    exit 1;
fi


PAIRED_OR_SINGLE_END="${1}";
BAM_LIST_FILENAME="${2}";
VCF_INPUT_FILENAME="${3}";
VCF_OUTPUT_FILENAME="${4}";
ASSAY_TYPE="${5}"
SINGLE_END="";

if [ "${PAIRED_OR_SINGLE_END}" = "single_end" ] ; then
    SINGLE_END="-singleEnd";
elif [ "${PAIRED_OR_SINGLE_END}" = "paired_end" ] ; then
    SINGLE_END="";
else
    printf 'Error: Specify read type: "single_end" or "paired_end"\n' > /dev/stderr;
    exit 1;
fi

MIN_INSERT=1
MAX_INSERT=270000000
if [ "${ASSAY_TYPE}" = "atac" ] ; then
    printf 'Assay type : ATAC-seq\n' > /dev/stderr;
    MIN_INSERT=37;
    MAX_INSERT=10000;
elif [ "${ASSAY_TYPE}" = "rna" ] ; then
    printf 'Assay type : RNA-seq\n' > /dev/stderr;
    MIN_INSERT=37;
    MAX_INSERT=2473537
fi

if [ ! -f "${BAM_LIST_FILENAME}" ] ; then
    printf 'Error: BAM list file "%s" does not exist.\n' "${BAM_LIST_FILENAME}" > /dev/stderr;
    exit 1;
fi

if [ ! -f "${VCF_INPUT_FILENAME}" ] ; then
    printf 'Error: VCF file "%s" does not exist.\n' "${VCF_INPUT_FILENAME}" > /dev/stderr;
    exit 1;
fi

if [ "${VCF_INPUT_FILENAME}" = "${VCF_INPUT_FILENAME%.vcf.gz}" ] ; then
    printf 'Error: VCF file "%s" must end with ".vcf.gz" (created with bgzip).\n' "${VCF_INPUT_FILENAME}" > /dev/stderr;
    exit 1;
fi


if [ -f "${VCF_OUTPUT_FILENAME}" ] ; then
    printf 'Error: VCF output file "%s" does already exist.\n' "${VCF_OUTPUT_FILENAME}" > /dev/stderr;
    exit 1;
fi


# Set ASVCF specific temp directory if the default temp dir is not big enough.
ASVCF_TMP_DIR="";
# Create directory for temporary files.
#   - if ${ASVCF_TMP_DIR} variable is not empty, create a subdir in that directory.
#   - else create a subdirectory in ${TMPDIR}.
#   - else if ${TMPDIR} is not set, create a subdirectory in /tmp.
ASVCF_CURRENT_TMP_DIR=$(mktemp -p "${ASVCF_TMP_DIR:=${TMPDIR}}" -d "ASVCF_${VCF_FILENAME##*/}_XXXXXX");

if [ ! -d "${ASVCF_CURRENT_TMP_DIR}" ] ; then
    printf 'Error: Unable to create temp directory "%s".\n' "${ASVCF_CURRENT_TMP_DIR}" > /dev/stderr;
    exit 1;
fi


ALL_ASC_FILENAME="${ASVCF_CURRENT_TMP_DIR}/all_asc.as.gz";


# Indexed array to keep track of the order of the BAM filenames as specified in ${BAM_LIST_FILENAME}.
declare -a BAM_FILES
# Indexed array to keep track of the order of the AS filenames.
declare -a AS_FILES

for BAM_FILENAME in $(cat "${BAM_LIST_FILENAME}") ; do
    if [ ! -f "${BAM_FILENAME}" ] ; then
        printf 'Error: BAM file "%s" does not exist.\n' "${BAM_FILENAME}" > /dev/stderr;
        exit 1;
    fi

    # Get the sample name by taking the basename of the BAM file withotu the ".bam" extension.
    SAMPLE_NAME="${BAM_FILENAME##*/}";
    SAMPLE_NAME="${SAMPLE_NAME%.bam}";

    AS_FILENAME="${ASVCF_CURRENT_TMP_DIR}/${SAMPLE_NAME}.as.gz";

    BAM_FILES+=("${BAM_FILENAME}");
    AS_FILES+=("${AS_FILENAME}");
done


# Count AS reads.
for CHR in $(tabix -l "${VCF_INPUT_FILENAME}") ; do
    # Create a temporary VCF file for the current chromosome only.
    TMP_VCF_FILENAME="${ASVCF_CURRENT_TMP_DIR}/${VCF_INPUT_FILENAME##*/}"

    tabix "${VCF_INPUT_FILENAME}" "${CHR}" \
        | cut -f 1-9 \
        | gzip \
        > "${TMP_VCF_FILENAME}";

    START=$(zcat "${TMP_VCF_FILENAME}" | head -n 1 | cut -f 2);
    END=$(zcat "${TMP_VCF_FILENAME}" | tail -n 1 | cut -f 2);

    echo "${CHR}:${START}-${END}";

    # Loop over each BAM file.
    for BAM_FILENAME_IDX in "${!BAM_FILES[@]}" ; do
        BAM_FILENAME="${BAM_FILES[${BAM_FILENAME_IDX}]}";
        AS_FILENAME="${AS_FILES[${BAM_FILENAME_IDX}]}";

        # Filter out secondary alignments from the BAM file.
        samtools view -F 0x0100 "${BAM_FILENAME}" "${CHR}:${START}-${END}" \
            | awk '$7=="="{print}' \
            | sort -k 1 \
            | "${RASQUALDIR}/src/ASVCF/qcFilterBam" stdin -skipMissing F -maxMismatch 3 -maxGapOpen 0 -maxBestHit 1 -minQual 10 -minInsert "${MIN_INSERT}" -maxInsert "${MAX_INSERT}" "${SINGLE_END}" \
            | sort -k 2 -n \
            | "${RASQUALDIR}/src/ASVCF/countAS" "${TMP_VCF_FILENAME}" \
            | awk -F '\t' '{ print $5 "," $6; }' \
            | gzip \
            >> "${AS_FILENAME}";
    done

    # Remove temporary VCF file.
    rm "${TMP_VCF_FILENAME}";
done


# Paste all AS counts in one file.
"${RASQUALDIR}/src/ASVCF/zpaste" "${AS_FILES[@]}" > "${ALL_ASC_FILENAME}";


# Remove all temporary AS files.
rm "${AS_FILES[@]}";


# merge AS counts and VCF.
"${RASQUALDIR}/src/ASVCF/pasteFiles" "${VCF_INPUT_FILENAME}" "${ALL_ASC_FILENAME}" \
    | bgzip \
    > "${VCF_OUTPUT_FILENAME}";


# Remove temporary files and directories.
rm "${ALL_ASC_FILENAME}";
rmdir "${ASVCF_CURRENT_TMP_DIR}";
