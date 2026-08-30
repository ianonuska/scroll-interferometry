# Reproducible environment for the winding meter (scroll-interferometry).
# Build:  docker build -t winding-meter .
# Demo:   docker run --rm -v $PWD/out:/out winding-meter \
#           python examples/run_demo.py            # streams PHerc1667 from S3
# Solve a custom slice: mount data and use winding_phase.winding_coordinate —
# see README "Run it" for the API.
FROM python:3.12-slim
RUN pip install --no-cache-dir \
    numpy scipy zarr s3fs matplotlib pillow tifffile pyamg
WORKDIR /app
COPY . /app
CMD ["python", "examples/run_demo.py"]
