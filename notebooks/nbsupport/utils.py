import json
def ppo(o):
    """pretty print object as json"""
    print(json.dumps(o.as_dict(), indent=2))
