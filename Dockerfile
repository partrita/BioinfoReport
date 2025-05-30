FROM ubuntu:24.04

# 1. 시스템 패키지 설치 (python, pip, LaTeX, latexmk, make 등)
RUN apt update && \
    apt install -y --no-install-recommends \
        python3-full \
        texlive-full \
        latexmk \
        make \
        git && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*


# 2. 작업 디렉토리 설정
WORKDIR /workspace

# 3. requirements.txt 복사 및 Python 의존성 설치
COPY requirements.txt .
RUN python3 -m venv /venv && \
    /venv/bin/pip3 install --upgrade pip && \
    /venv/bin/pip3 install --no-cache-dir -r requirements.txt

# 4. 전체 프로젝트 복사
COPY . .

# Ensure venv is activated for subsequent commands
ENV PATH="/venv/bin:$PATH"

# 5. 기본 명령어를 bash로 설정
CMD ["/bin/bash"]
