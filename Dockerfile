#Dockerfile

FROM python:3.10.1

ENV USER=streamlit
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/app

# Update system and install dependencies.
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    software-properties-common

RUN apt-get install libxrender1
COPY . $HOME/app/
RUN pip install --no-cache-dir -r requirements.txt

RUN ["chmod", "+x", "start-script.sh"]

EXPOSE 8501

CMD ./start-script.sh