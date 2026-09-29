FROM public.ecr.aws/d3j8x8q7/olympus-base-python:latest

WORKDIR /app
COPY . .

RUN pip install -e .
RUN pip install rdflib==7.6.0 reportlab==5.0.1 networkx==3.7 scipy==1.18.1 igraph==1.0.0 matplotlib==3.11.2 mmtf-python==1.1.3
CMD ["/bin/bash"]
