# Используем базовый образ с tidyverse
FROM rocker/tidyverse:latest

# --- Системные зависимости ---
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk40 \
    && rm -rf /var/lib/apt/lists/*

# --- Устанавливаем CRAN пакеты ---
RUN R -e "install.packages(c('dplyr', 'stringr', 'readr', 'ggplot2', 'optparse'), repos='https://cloud.r-project.org')"

# --- Устанавливаем Bioconductor пакеты ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('VariantAnnotation'), update=FALSE, ask=FALSE)"

# --- Копируем R-скрипты (VCF статистика) ---
COPY scripts/ /usr/local/bin/

# --- Делаем скрипты исполняемыми ---
RUN chmod +x /usr/local/bin/*.R

# --- Создаем алиасы без .R для удобного вызова ---
RUN for f in /usr/local/bin/*.R; do ln -s "$f" "${f%.R}"; done

# --- Рабочая директория ---
WORKDIR /usr/local/bin

# --- По умолчанию интерактивный shell ---
ENTRYPOINT ["/bin/bash"]