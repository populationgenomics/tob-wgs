FROM australia-southeast1-docker.pkg.dev/peter-dev-302805/test/hailr:hail0.2.73_r4.1
LABEL maintainer="https://github.com/pdiakumis"

# Concordance scripts + Rmd
ENV SCRIPT_PATH "/home/scripts"
ENV PATH "${SCRIPT_PATH}:${PATH}"
WORKDIR ${SCRIPT_PATH}
COPY concordance concordance.Rmd ${SCRIPT_PATH}/
ENV TZ "Australia/Melbourne"

ENTRYPOINT [ "concordance" ]
