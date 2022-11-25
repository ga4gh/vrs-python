FROM python:3.10
COPY . .
RUN pip install -e .[dev,extras]
RUN mkdir /data
CMD ["bash"]
