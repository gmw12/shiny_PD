spd-app
=======

OpenShift deployment of spd application using shiny server
## Contents

- Storage.yaml: Configuration for requesting persistent storage to house data
- Build.yaml: Build configurations and image streams for base image and application image
- Deployment.yaml: Deployment configuration to run application pods

We have a two step build process(base image and application image) so code changes can be deployed quickly.
Otherwise code changes would require waiting for all requirements to be installed.

## Deployment

1. Create shared storage `oc create -f Storage.yaml`
3. Create Build and image configurations `oc create -f Build.yaml`
4. Deploy application: `oc create -f Deployment.yaml`
5. Ask openshift for the application route: `oc get route spd-shiny-route`
