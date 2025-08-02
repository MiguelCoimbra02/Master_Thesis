FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml first to leverage Docker cache
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml

# Ensure conda environment is activated by default for RUN commands
SHELL ["conda", "run", "-n", "viz_tool", "/bin/bash", "-c"]

# Copy application code
COPY . .

# Expose port 5001 for Flask app
EXPOSE 5001

# Use conda run to activate environment and run app.py by default
CMD ["conda", "run", "--no-capture-output", "-n", "viz_tool", "python", "-u", "app.py"]