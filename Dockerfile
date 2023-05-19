FROM base-shiny:latest
ADD . /code
RUN chmod 777 /code/restart.txt
