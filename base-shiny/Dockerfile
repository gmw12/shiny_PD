FROM rocker/shiny:4.1.1 #-verse:4.0.2

RUN apt-get update
RUN apt-get install build-essential libgl1-mesa-dev libglu1-mesa-dev libv8-dev libglpk-dev libbz2-dev liblzma-dev -y

ADD ./Shiny_install.R /code/Shiny_install.R
RUN Rscript /code/Shiny_install.R

COPY base-shiny/shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN chown -R shiny /var/lib/shiny-server/

# OpenShift gives a random uid for the user and some programs try to find a username from the /etc/passwd.
# Let user to fix it, but obviously this shouldn't be run outside OpenShift
RUN chmod ug+rw /etc/passwd
COPY base-shiny/fix-username.sh /fix-username.sh
COPY base-shiny/shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod a+rx /usr/bin/shiny-server.sh

# Make sure the directory for individual app logs exists and is usable
RUN chmod -R a+rwX /var/log/shiny-server
RUN chmod -R a+rwX /var/lib/shiny-server

CMD /usr/bin/shiny-server.sh
