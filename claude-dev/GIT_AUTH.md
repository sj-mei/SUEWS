# Git Authentication in Claude Sandbox

## Recommended: GitHub CLI (Simple & Secure)

1. **Inside the container**, run:
   ```bash
   gh auth login
   ```

2. Choose:
   - GitHub.com
   - HTTPS
   - Login with browser

3. You're authenticated! Git operations will now work.

**Benefits:**
- ✅ No SSH keys in container
- ✅ Works on all platforms
- ✅ Tokens can be revoked
- ✅ Fine-grained permissions
- ✅ Perfect for team use

## Alternative: SSH Agent Forwarding (Advanced)

### Linux
```bash
# Add to your shell profile
./claude.sh start my-workspace \
  -v $SSH_AUTH_SOCK:$SSH_AUTH_SOCK \
  -e SSH_AUTH_SOCK=$SSH_AUTH_SOCK
```

### macOS
```bash
# Docker Desktop provides a special socket
./claude.sh start my-workspace \
  -v /run/host-services/ssh-auth.sock:/run/host-services/ssh-auth.sock \
  -e SSH_AUTH_SOCK=/run/host-services/ssh-auth.sock
```

### Windows (WSL2)
```bash
# In WSL2, use Linux method
./claude.sh start my-workspace \
  -v $SSH_AUTH_SOCK:$SSH_AUTH_SOCK \
  -e SSH_AUTH_SOCK=$SSH_AUTH_SOCK
```

## Why NOT Mount SSH Keys?

Mounting SSH keys directly (`~/.ssh:/tmp/.ssh`) is:
- ❌ Security risk - keys persist in container
- ❌ Keys visible in `docker inspect`
- ❌ Difficult to revoke access
- ❌ Against Docker security best practices

## For CI/CD

Use GitHub Actions secrets:
```yaml
- name: Authenticate
  run: |
    echo "${{ secrets.GITHUB_TOKEN }}" | gh auth login --with-token
```

## Troubleshooting

If `gh auth login` fails:
1. Ensure you have internet access
2. Try `gh auth logout` first
3. Check `gh auth status`

If SSH agent forwarding fails:
1. Ensure ssh-agent is running: `ssh-add -l`
2. Add your key: `ssh-add ~/.ssh/id_ed25519`
3. Check socket exists: `echo $SSH_AUTH_SOCK`