FROM base-shiny:latest
ADD . /code
RUN chmod 777 /code/restart.txt
RUN chown shiny /code/restart.txt
