
# Development Analysis Report

[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://partrita.github.io)

이 프로젝트는 **단백질 의약품 개발** 과정에서 자주 사용되는 바이오인포매틱스 도구와 분석 파이프라인을 재현 가능하고 체계적으로 관리하기 위해 시작되었습니다.

실제 현장에서는 다양한 도구들이 경험과 직관에 의존해 사용되는 경우가 많아, 표준화된 분석 파이프라인의 필요성을 느꼈습니다. [Amber Biology](https://www.amberbiology.com/)의 항체 분석 데모 리포트를 보고, Jupyter Notebook 기반으로 누구나 재현 가능한 분석 환경을 만들고자 본 프로젝트를 시작하였습니다.

## 주요 특징

- Jupyter Book 기반의 문서화 및 PDF 리포트 자동 생성
- Docker를 활용한 일관된 개발/분석 환경 제공
- 바이오인포매틱스 및 수치해석 예제 코드 수록


## 환경 구성 및 사용법

### 1. Docker 이미지 빌드

```bash
docker build -t jb-latex-env .
```

### 2. 컨테이너 실행 (현재 디렉토리 마운트)

```bash
docker run --rm -it --user $(id -u):$(id -g) -v $(pwd):/workspace jb-latex-env
```

### 3. 컨테이너 내에서 Jupyter Book 빌드

```bash
jupyter-book build demo-report/Aflibercept --all
```

#### PDF 파일 생성 및 이동

```bash
jupyter-book build demo-report/Aflibercept --builder pdflatex
mv demo-report/Aflibercept/_build/latex/*.pdf demo-report/Aflibercept/
```

## 라이선스

이 프로젝트는 MIT 라이선스를 따릅니다.  
자세한 내용은 [LICENSE](LICENSE) 파일을 참고하세요.

---

## 참고

- [Jupyter Book 공식 문서](https://jupyterbook.org)
- [Amber Biology](https://www.amberbiology.com)
- [partrita.github.io](https://partrita.github.io)