from torch import device, cuda

def get_device(dev: str | device | None, verbose: bool = False) -> device:
    """ (from Jeff's device_handling: 20251003)
    Resolve user's requested device.
    Args:
        dev (str | device | None): Requested device. If a torch.device, it will be returned
            unchanged. If None, will default to CUDA:0 if CUDA is available, or else will use
            cpu. String values are any recognized by torch, i.e. cpu, cuda with or without
            ordinal, or other device types. We do not attempt to ensure the requested
            device is actually available on the system.
        verbose (bool, optional): If True, function will describe the request and response.
            Defaults to False.
    Returns:
        device: A torch.device initialized as requested.
    """
    if isinstance(dev, device):
        if verbose:
            if dev.type == 'cuda':
                print(f"Using prepared device {cuda.get_device_properties(dev)}")
            else:
                print(f"Using prepared device of type {dev.type}")
        return dev
    elif dev is None:
        if cuda.is_available():
            if verbose:
                print(f"No device requested, defaulting to available cuda.")
            return device('cuda')
        if verbose:
            print(f"No device requested, CUDA unavailable, defaulting to CPU")
        return device('cpu')
    if not cuda.is_available():
        if verbose:
            print(f"Device {dev} was requested but CUDA is not available, using CPU.")
        return device('cpu')
    if verbose:
        print(f"Using requested device {dev}")
    return device(dev)

