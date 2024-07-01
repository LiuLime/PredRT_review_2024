"""
Send Http requirements
@ Liu
"""

import socket
import time
import urllib.error

import urllib3
from urllib3.exceptions import InsecureRequestWarning
import requests
from urllib.error import URLError, HTTPError


class ErrorHandle:
    def __init__(self):
        pass

    def error_handle(self, e, entity, attempt):
        match e:
            case urllib.error.URLError() as url_err:
                if isinstance(url_err.reason, socket.timeout):
                    print(entity, "time out")
                    return self.retry(attempt)

            case _:
                print(entity, "unknown error")
                return False

    def retry(self, attempt, max_retries=5, retry_delay=3):
        if attempt < max_retries:
            print(f"{attempt}/{max_retries} retry....")
            time.sleep(retry_delay)
            return True
        else:
            print(f"failed {max_retries} times, break process")
            return False


def send_url(url, headers=None, verify=False, max_retries=2, rate_limit=0.5) -> dict | None | str:
    """
    send http requests
    :param url:
    :param headers:
    :param verify:
    :param max_retries:
    :return:
    """
    urllib3.disable_warnings(InsecureRequestWarning)
    time.sleep(rate_limit)
    entry = None
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url, headers=headers, verify=verify)
            response.raise_for_status()  # This will raise HTTPError for bad responses
            print(f"{url} response")
            entry = response.text
            break  # Success, exit the loop
        except HTTPError as e:
            print(f"{url} HTTP error: {e}")
            break  # Break on HTTP errors
        except ConnectionError as e:
            retries += 1
            print(f"{url} connection failed: {e}, attempt {retries} times")
            time.sleep(3 ** retries)  # Exponential backoff
        except Exception as e:
            print(f"{url} failed for {e}")
            break
    return entry