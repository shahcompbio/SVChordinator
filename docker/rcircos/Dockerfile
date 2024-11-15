# Use tidyverse base
FROM rocker/tidyverse:4.4.2
LABEL authors="asherpreskasteinberg"

# install r-circos
RUN R -e 'install.packages("RCircos")'
# Set default to bash for script execution flexibility
CMD ["bash"]

