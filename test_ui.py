import requests
r = requests.get('http://127.0.0.1:5000/')
print('Status:', r.status_code, 'Size:', len(r.text))
print('Has progress-bar:', 'progress-bar-wrapper' in r.text)
print('Has help-overlay:', 'help-overlay' in r.text)
print('Has help-btn:', 'help-btn' in r.text)
print('Has old loading div:', 'id="loading"' in r.text)
