# Use Python base image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install Python dependencies
COPY app/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app code
COPY app /app
COPY frontend /app/frontend

# Expose port for uvicorn
EXPOSE 8000

# Run FastAPI
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
