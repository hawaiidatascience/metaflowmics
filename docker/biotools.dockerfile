FROM quay.io/biocontainers/seqtk:1.3--h5bf99c6_3
FROM biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1

COPY --from=0 /usr/local/bin/seqtk /usr/bin/
