# Используем базовый образ с tidyverse
FROM rocker/tidyverse:latest

# --- Системные зависимости ---
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk40 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# --- Устанавливаем R пакеты (CRAN) ---
RUN R -e "install.packages(c('dplyr', 'stringr', 'readr', 'ggplot2', 'optparse'), repos='https://cloud.r-project.org')"

# --- Устанавливаем Bioconductor пакеты ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('VariantAnnotation'), update=FALSE, ask=FALSE)"

# --- Устанавливаем Python пакеты для скрипта ---
RUN pip3 install cyvcf2 pandas

# --- Копируем скрипты в контейнер ---
COPY scripts/ /usr/local/bin/

# --- Делаем скрипты исполняемыми ---
RUN chmod +x /usr/local/bin/*.R
RUN chmod +x /usr/local/bin/vcf_to_gt_table.py

# --- Создаём алиас без расширения для Python скрипта ---
RUN ln -s /usr/local/bin/vcf_to_gt_table.py /usr/local/bin/vcf_to_gt_table

# --- Рабочая директория ---
WORKDIR /usr/local/bin

# --- По умолчанию интерактивный shell ---
ENTRYPOINT ["/bin/bash"]