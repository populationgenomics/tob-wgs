FROM pdiakumis/hailr:r4.1_hail0.2.70
LABEL maintainer="https://github.com/pdiakumis"

# Concordance scripts + Rmd
ENV SCRIPT_PATH "/home/scripts"
ENV PATH "${SCRIPT_PATH}:${PATH}"
WORKDIR ${SCRIPT_PATH}
COPY concordance concordance.Rmd .
ENV TZ "Australia/Melbourne"

ENTRYPOINT [ "concordance" ]
