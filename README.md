# georges
Georges' the lemur opinionated particle accelerator modeling Python package.

<img src="https://github.com/chernals/georges/blob/master/georges.png" width="300" />

## Docker image
A Docker image is made available to provide an easy access to a complete Jupyter Notebook + madx + georges environment.
 
 Use  the *Dockerfile* to build the image:
 
 ```
 docker build .
```

or, to register the image as well:

```
docker build -t username/georges .
```

You can run a container with

```
docker run -it username/georges
```

then connect to [http://localhost:8888](http://localhost:8888 "Jupyter Notebook") to access the Jupyter Notebook interface.

The image includes a complete Anaconda Python3 environment with the most common packages. 
The latest *MAD-X* development release is available in */usr/local/bin/madx*.