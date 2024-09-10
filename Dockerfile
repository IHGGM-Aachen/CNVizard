# Use an official Python runtime as a parent image
FROM python:3.12.4-slim

# Install non-python package dependencies
RUN apt-get update && apt-get install -y \
    tabix \
    && apt-get clean

# Set the working directory in the container
WORKDIR /app

# Copy only the necessary files for installing the application
COPY setup.py MANIFEST.in default.env README.md ./
COPY cnvizard ./cnvizard

# Install the package using setup.py
RUN pip install --no-cache-dir .

# Expose the port the app runs on
EXPOSE 8501

# Run Streamlit
CMD ["streamlit", "run", "cnvizard/app.py"]
