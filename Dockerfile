FROM python:3.12-slim

LABEL org.opencontainers.image.title="BOTAS"
LABEL org.opencontainers.image.description="Bacterial Operon-aware Transcriptome Alignment System"
LABEL org.opencontainers.image.source="https://github.com/clabe-wekesa/botas"
LABEL org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /opt/botas

COPY . .

RUN python -m pip install --upgrade pip \
    && python -m pip install ".[progress]"

WORKDIR /work

ENTRYPOINT ["botas"]
CMD ["--help"]
