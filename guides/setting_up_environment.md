# Setting Up the BulkAlphaFold3 Environment

A guide to setting up the environment for BulkAlphaFold3, tailored to the LambdaLabs environment. This guide can also be used to just set up AlphaFold3 as well. This is assuming you have a LambdaLabs account with billing properly set up and are ready to launch an instance.

**In total, this guide should take no longer than an hour to complete, assuming you have model weights for AlphaFold3 and a solid internet connection. A bulk of this guide is just properly setting up the instance and installing AlphaFold3.**

## Setting Up An Instance with AlphaFold3

After setting up a LambdaLabs account, launch a new instance. In this guide, we will use 1x GH200 instance with 96GB of GPU memory, sufficient for running predictions of up to length [XXXX] residues.

From a terminal in the instance (which can be accessed via SSH or just by launching their convienient JupyterLab web interface), follow the following steps.

### Install Docker and Docker GPU Support

Follow the `Installing Docker` instructions on the AlphaFold3 [installation page](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#installing-docker). This includes `Installing Docker on Host` and `Enabling Rootless Docker`.

For reference, here are the commands you will need to run:

```bash
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo docker run hello-world

sudo apt-get install -y uidmap systemd-container
sudo machinectl shell $(whoami)@ /bin/bash -c 'dockerd-rootless-setuptool.sh install && sudo loginctl enable-linger $(whoami) && DOCKER_HOST=unix:///run/user/1001/docker.sock docker context use rootless'
```

**After, you do not need to go through the complete `Installing GPU Support` section, as the LambdaLabs instance already has GPU support enabled.** Instead, you only need to do the `Installing NVIDIA Support for Docker` section, and you should not need to run `systemctl --user restart docker`, i.e. run: 

```bash
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
nvidia-ctk runtime configure --runtime=docker --config=$HOME/.config/docker/daemon.json
# systemctl --user restart docker # Will not work on LambdaLabs instance, skipping is just fine
sudo nvidia-ctk config --set nvidia-container-cli.no-cgroups --in-place
```

The following command should return the NVIDIA driver version and GPU information, confirming that Docker is set up correctly:

```
docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu22.04 nvidia-smi
```

### Installing AlphaFold3

Run the following commands to install AlphaFold3 from Github:

```bash
cd ~/
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3/
bash fetch_databases.sh
```

The following command may take awhile, and so in the meantime, you can open a new terminal or tab in the JupyterLab interface to continue with building the Docker image. 

### Building the Docker Image

Because the LambdaLabs instance has an ARM64 architecture, we need to slightly modify the Dockerfile to build the image correctly (see https://github.com/google-deepmind/alphafold3/issues/422 for details).

If you are using x86_64 architecture, **you can skip this step** and just run the `docker build` command below.

```bash
wget -O /home/ubuntu/alphafold3/docker/Dockerfile https://raw.githubusercontent.com/daniel-ji/alphafold3/refs/heads/main/docker/Dockerfile
```

Now, build the Docker image with the following command:

```bash
cd ~/alphafold3/
docker build -t alphafold3 -f docker/Dockerfile .
```


### Uploading AlphaFold3 Model Weights

After obtaining AlphaFold3 model weights (refer to https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#obtaining-model-parameters or contact Daniel at daji@ucsd.edu), upload them to the instance. It should be a single file titled `af3.bin.zst` (~1GB). You can use `scp` or any other file transfer method you prefer (or use the JupyterLab interface to upload the weights directly, but this may be unstable at times). **In these guides, we will assume the file is uploaded to the home directory of the instance.**

## Installing BulkAlphaFold3

Now that you have AlphaFold3 installed and the model weights uploaded, you can install BulkAlphaFold3.

Simply clone the repository and install the required Python packages:

```bash
cd ~/
git clone git@github.com:daniel-ji/bulkalphafold3.git
cd bulkalphafold3/
pip install -r requirements.txt
```

And now you are ready to run BulkAlphaFold3. Installation is complete!