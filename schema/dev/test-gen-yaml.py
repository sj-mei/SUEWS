# from jsf import JSF
import json
from hypothesis import given
from hypothesis_jsonschema import from_schema

# faker = JSF(
#     {
#         "type": "object",
#         "properties": {
#             "name": {"type": "string", "$provider": "faker.name"},
#             "email": {"type": "string", "$provider": "faker.email"},
#         },
#         "required": ["name", "email"],
#     }
# )
with open("suews-var-def.json") as f:
    user_schema = json.load(f)

# print(user_schema)


@given(from_schema(user_schema))
def test_integers(value):
    print(len(value))

    assert isinstance(value[2]['price'], int)

# faker = JSF(user_schema)

# test_integers()


if __name__ == "__main__":
    test_integers()